//
//  enquserver-node.cpp
//  NS3-cc
//
//  Created by 王钧 on 2023/9/26.
//

#include "ns3/ipv4.h"
#include "ns3/packet.h"
#include "ns3/ipv4-header.h"
#include "ns3/pause-header.h"
#include "ns3/flow-id-tag.h"
#include "ns3/boolean.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"
#include "enquserver-node.h"
#include "qbb-net-device.h"
#include "ppp-header.h"
#include "ns3/int-header.h"
#include <cmath>

namespace ns3 {

TypeId EnquserverNode::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::EnquserverNode")
    .SetParent<Node> ()
    .AddConstructor<EnquserverNode> ()
    .AddAttribute("EcnEnabled",
            "Enable ECN marking.",
            BooleanValue(false),
            MakeBooleanAccessor(&EnquserverNode::m_ecnEnabled),
            MakeBooleanChecker())
    .AddAttribute("CcMode",
            "CC mode.",
            UintegerValue(0),
            MakeUintegerAccessor(&EnquserverNode::m_ccMode),
            MakeUintegerChecker<uint32_t>())
    .AddAttribute("AckHighPrio",
            "Set high priority for ACK/NACK or not",
            UintegerValue(0),
            MakeUintegerAccessor(&EnquserverNode::m_ackHighPrio),
            MakeUintegerChecker<uint32_t>())
    .AddAttribute("MaxRtt",
            "Max Rtt of the network",
            UintegerValue(9000),
            MakeUintegerAccessor(&EnquserverNode::m_maxRtt),
            MakeUintegerChecker<uint32_t>())
  ;
  return tid;
}

EnquserverNode::EnquserverNode(){
    m_ecmpSeed = m_id;
    m_node_type = 1;
    m_mmu = CreateObject<SwitchMmu>();
    for (uint32_t i = 0; i < pCnt; i++)
        for (uint32_t j = 0; j < pCnt; j++)
            for (uint32_t k = 0; k < qCnt; k++)
                m_bytes[i][j][k] = 0;
    for (uint32_t i = 0; i < pCnt; i++)
        m_txBytes[i] = 0;
    for (uint32_t i = 0; i < pCnt; i++)
        m_lastPktSize[i] = m_lastPktTs[i] = 0;
    for (uint32_t i = 0; i < pCnt; i++)
        m_u[i] = 0;
}

int EnquserverNode::GetOutDev(Ptr<const Packet>p, CustomHeader &ch){
    // look up entries
    auto entry = m_routerMap.find(ch.dip);

    // no matching entry
    if (entry == m_routerMap.end())
        return -1;

    // entry found
    auto idx = entry->second;


    return idx;
}

void EnquserverNode::CheckAndSendPfc(uint32_t inDev, uint32_t qIndex){
    Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
    if (m_mmu->CheckShouldPause(inDev, qIndex)){
        device->SendPfc(qIndex, 0);
        m_mmu->SetPause(inDev, qIndex);
    }
}
void EnquserverNode::CheckAndSendResume(uint32_t inDev, uint32_t qIndex){
    Ptr<QbbNetDevice> device = DynamicCast<QbbNetDevice>(m_devices[inDev]);
    if (m_mmu->CheckShouldResume(inDev, qIndex)){
        device->SendPfc(qIndex, 1);
        m_mmu->SetResume(inDev, qIndex);
    }
}

void EnquserverNode::SendToDev(Ptr<Packet>p, CustomHeader &ch){
    int idx = GetOutDev(p, ch);
    if (idx >= 0){
        NS_ASSERT_MSG(m_devices[idx]->IsLinkUp(), "The routing table look up should return link that is up");

        // determine the qIndex
        uint32_t qIndex = 0;//出入口服务器接收到的一定是ACK包，设置优先级为0

        m_devices[idx]->SwitchSend(qIndex, p);
    }else
        return; // Drop
}




//生成共享链路表操作

void EnquserverNode::GetShareTable(Ptr<const Packet>p, CustomHeader &ch){
    if (ch.tcp.fin==1){//接收到fin标识位为1，说明该流结束，此时将该数据包中的{sip,dip,sport,dport}对应的共享链路表中的表项全部删除
        for (const m_sharedTableEntry& p : m_sharedTable) {
            p.flowInfos.erase(std::remove_if(p.flowInfos.begin(), p.flowInfos.end(),
                                            [ch.sip,ch.dip,ch.tcp.sport,ch.tcp.dport](const flowInfo& f) {
                                                return f.sip == ch.dip && f.dip == ch.sip&& f.dport == ch.tcp.sport&& f.sport == ch.tcp.dport;
                                            }),
                              p.flowInfos.end());
            }
    }
    else{//获取接收到的ack包中的路由id和port信息，在共享链路表对应的表项中查找，若没有，则直接添加
        bool found = false;
           for (m_sharedTableEntry& p : m_sharedTable) {
               if (p.rid == ch.ih.iinfo.id && p.port == ch.ih.iinfo.port) {
                   found = true;
                   p.flowInfos.push_back({ch.dip,ch.sip,ch.tcp.dport,ch.tcp.sport});//将该数据包的四元组信息添加到对应的表项中
                   break; // 如果找到了，跳出循环
               }
           }

           // 如果没找到，添加到向量中
           if (!found) {
               std::vector<flowInfo> info;
               info.push_back({ch.dip,ch.sip,ch.tcp.dport,ch.tcp.sport}); //将四元组信息添加到Info中
               m_sharedTable.push_back({ch.ih.iinfo.id,ch.ih.iinfo.port,info})；
           }
    }
}


//对携带链路信息的数据包中的信息和共享链路表进行查找匹配，返回HeaderLinkInfo结构体类型中的数据
void EnquserverNode::MatchSharedTableSendToRelatedSender(Ptr<NetDevice> device, Ptr<Packet>p, CustomHeader &ch){
    GetShareTable(packet, ch)；
    
    std::vector<m_sharedTableEntry> matchedEntries;//根据从数据包获取的路由器二元组信息，和共享链路表比配，获取数据包中路由节点二元组对应的所有主机地址四元组
    for (int i = 0; i < ch.tcp.ih.hinfo.depthNum; ++i) {
        for (const auto& sharedEntry : m_sharedTable) {
            if (ch.tcp.ih.dinfo[i].iinfo.id == sharedEntry.rid && ch.tcp.ih.dinfo[i].iinfo.port == sharedEntry.port) {
                matchedEntries.push_back(sharedEntry);
            }
        }
    }
    for (int i = 0; i < ch.tcp.ih.hinfo.ratioNum; ++i) {
        for (const auto& sharedEntry : m_sharedTable) {
            if (ch.tcp.ih.rinfo[i].iinfo.id == sharedEntry.rid && ch.tcp.ih.rinfo[i].iinfo.port == sharedEntry.port) {
                matchedEntries.push_back(sharedEntry);
            }
        }
    }
    
    std::vector<relatedSenderHeaderInfo> relatedSenderHeaderInfos;
    for (const auto& sharedEntry : matchedEntries) {
        for (const auto& flow : sharedEntry.flowInfos) {
            relatedSenderHeaderInfo relatedInfo;

            // 设置 fInfo
            relatedInfo.fInfo = flow;

            // 查找是否已经存在对应的 rIdAndPort
            auto it = std::find_if(
                relatedSenderHeaderInfos.begin(),
                relatedSenderHeaderInfos.end(),
                [flow](const relatedSenderHeaderInfo& info) {
                    return info.fInfo == flow;
                }
            );

            // 如果存在，直接添加 rIdAndPort，否则添加新的记录
            if (it != relatedSenderHeaderInfos.end()) {
                
                it->rIdAndPort.push_back({sharedEntry.rid, sharedEntry.port});
            } else {
                relatedInfo.rIdAndPort.push_back({sharedEntry.rid, sharedEntry.port});
                relatedSenderHeaderInfos.push_back(relatedInfo);
            }
        }
    }
    for (const auto& info : relatedSenderHeaderInfos) { //需要加判断，如果sip。。。。==原数据包中的sip。。。。，则直接转发
        if (info.fInfo.sip == ch.dip && info.fInfo.dip == ch.sip && info.fInfo.sport == ch.tcp.dport && info.fInfo.dport == ch.tcp.sport) {
            SendToDev(packet,info.fInfo.sip); //直接转发
        }else{
            encHeader encH;
    //        seqh.SetSeq(rxQp->ReceiverNextExpectedSeq);
            encH.SetFlags(1);
    //        seqh.SetPG(ch.udp.pg);
            encH.SetSport(info.fInfo.dport);
            encH.SetDport(info.fInfo.sport);
            
            MyIntHeader ih;
            for (const auto& pair : info.rIdAndPort) {
                bool isFind = 0;
                for (int i = 0; i < ch.tcp.ih.hinfo.depthNum; ++i){
                    if (pair.first == ch.tcp.ih.dinfo[i].iinfo.id && pair.second == ch.tcp.ih.dinfo[i].iinfo.port) {
                        ih.PushDepth(ch.tcp.ih.dinfo[i].iinfo.id,ch.tcp.ih.dinfo[i].iinfo.port,ch.tcp.ih.dinfo[i].depth,ch.tcp.ih.dinfo[i].ts,ch.tcp.ih.dinfo[i].maxRate);
                        isFind = 1;
                        break;
                    }
                }
                if (isFind == 1) {
                    break;
                }
                for (int i = 0; i < ch.tcp.ih.hinfo.ratioNum; ++i){
                    if (pair.first == ch.tcp.ih.rinfo[i].iinfo.id && pair.second == ch.tcp.ih.rinfo[i].iinfo.port) {
                        ih.PushDepth(ch.tcp.ih.rinfo[i].iinfo.id,ch.tcp.ih.rinfo[i].iinfo.port,ch.tcp.ih.rinfo[i].ratio,ch.tcp.ih.rinfo[i].ts,ch.tcp.ih.rinfo[i].maxRate);
                        break;
                    }
                }
            }
            encH.SetMyIntHeader(ih);
    //        if (ecnbits)
    //            seqh.SetCnp();

            Ptr<Packet> newp = Create<Packet>(0);
            newp->AddHeader(encH); //将ppp头部的上述信息写入到buffer中，方便后续在receive数据包时，ch从buffer中读取

            Ipv4Header head;    // Prepare IPv4 header
            head.SetDestination(Ipv4Address(info.fInfo.sip));
            head.SetSource(Ipv4Address(info.fInfo.dip));
            head.SetProtocol(0xFC); //ack=0xFC nack=0xFD
            head.SetTtl(64);
            head.SetPayloadSize(newp->GetSize());
    //        head.SetIdentification(rxQp->m_ipid++);

            newp->AddHeader(head);
            PppHeader ppp;
            ppp.SetProtocol (0x0021);//IPv4
            newp->AddHeader (ppp);
            CustomHeader ch(CustomHeader::L2_Header | CustomHeader::L3_Header | CustomHeader::L4_Header);
            p->PeekHeader(ch);
    //        AddHeader(newp, 0x800);    // Attach PPP header
            
//            CustomHeader ch(CustomHeader::L2_Header | CustomHeader::L3_Header | CustomHeader::L4_Header);
//            ch.getInt = 1; // parse INT header
//            newp->PeekHeader(ch); //把packet中的相关信息read到ch中
            SendToDev(newp, ch)
//            uint32_t nic_idx = GetNicIdxOfRxQp(info.fInfo.dip);
//            m_nic[nic_idx].dev->RdmaEnqueueHighPrioQ(newp);
//            m_nic[nic_idx].dev->TriggerTransmit();
        }
        
    }


}
//uint32_t EnquserverNode::GetNicIdxOfRxQp(uint32_t dip){
//    auto &v = m_rtTable[dip];
//    if (v.size() > 0){
//        return v[q->GetHash() % v.size()];
//    }else{
//        NS_ASSERT_MSG(false, "We assume at least one NIC is alive");
//    }
//}
//


uint32_t EnquserverNode::EcmpHash(const uint8_t* key, size_t len, uint32_t seed) {
  uint32_t h = seed;
  if (len > 3) {
    const uint32_t* key_x4 = (const uint32_t*) key;
    size_t i = len >> 2;
    do {
      uint32_t k = *key_x4++;
      k *= 0xcc9e2d51;
      k = (k << 15) | (k >> 17);
      k *= 0x1b873593;
      h ^= k;
      h = (h << 13) | (h >> 19);
      h += (h << 2) + 0xe6546b64;
    } while (--i);
    key = (const uint8_t*) key_x4;
  }
  if (len & 3) {
    size_t i = len & 3;
    uint32_t k = 0;
    key = &key[i - 1];
    do {
      k <<= 8;
      k |= *key--;
    } while (--i);
    k *= 0xcc9e2d51;
    k = (k << 15) | (k >> 17);
    k *= 0x1b873593;
    h ^= k;
  }
  h ^= len;
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

void EnquserverNode::SetEcmpSeed(uint32_t seed){
    m_ecmpSeed = seed;
}

void EnquserverNode::AddTableEntry(Ipv4Address &dstAddr, uint32_t intf_idx){
    uint32_t dip = dstAddr.Get();
    uint32_t outPort = intf_idx;
    m_routerMap[dip] = outPort;
}

void EnquserverNode::ClearTable(){
    m_routerMap.clear();
}

// This function can only be called in switch mode
//bool EnquserverNode::SwitchReceiveFromDevice(Ptr<NetDevice> device, Ptr<Packet> packet, CustomHeader &ch){
//    GetShareTable(packet, ch)；
//    SendToDev(packet, ch);
//    return true;
//}

//void EnquserverNode::EnquServerNotifyAllHost(uint32_t ifIndex, uint32_t qIndex, Ptr<Packet> p){ //移除头部ppp，移除IPV4头部，
//    FlowIdTag t;
//    p->PeekPacketTag(t);
//    if (qIndex != 0){
//        uint32_t inDev = t.GetFlowId();
//        m_mmu->RemoveFromIngressAdmission(inDev, qIndex, p->GetSize());
//        m_mmu->RemoveFromEgressAdmission(ifIndex, qIndex, p->GetSize());
//        m_bytes[inDev][ifIndex][qIndex] -= p->GetSize();
//        if (m_ecnEnabled){
//            bool egressCongested = m_mmu->ShouldSendCN(ifIndex, qIndex);
//            if (egressCongested){
//                PppHeader ppp;
//                Ipv4Header h;
//                p->RemoveHeader(ppp);
//                p->RemoveHeader(h);
//                h.SetEcn((Ipv4Header::EcnType)0x03);
//                p->AddHeader(h);
//                p->AddHeader(ppp);
//            }
//        }
//        //CheckAndSendPfc(inDev, qIndex);
//        CheckAndSendResume(inDev, qIndex);
//    }
//    if (1){
//        uint8_t* buf = p->GetBuffer();
//        if (buf[PppHeader::GetStaticSize() + 9] == 0x11){ // udp packet
//            IntHeader *ih = (IntHeader*)&buf[PppHeader::GetStaticSize() + 20 + 8 + 6]; // ppp, ip, udp, SeqTs, INT
//            Ptr<QbbNetDevice> dev = DynamicCast<QbbNetDevice>(m_devices[ifIndex]);
//            if (m_ccMode == 3){ // HPCC
//                ih->PushHop(Simulator::Now().GetTimeStep(), m_txBytes[ifIndex], dev->GetQueue()->GetNBytesTotal(), dev->GetDataRate().GetBitRate());
//            }else if (m_ccMode == 10){ // HPCC-PINT
//                uint64_t t = Simulator::Now().GetTimeStep();
//                uint64_t dt = t - m_lastPktTs[ifIndex];
//                if (dt > m_maxRtt)
//                    dt = m_maxRtt;
//                uint64_t B = dev->GetDataRate().GetBitRate() / 8; //Bps
//                uint64_t qlen = dev->GetQueue()->GetNBytesTotal();
//                double newU;
//
//                /**************************
//                 * approximate calc
//                 *************************/
//                int b = 20, m = 16, l = 20; // see log2apprx's paremeters
//                int sft = logres_shift(b,l);
//                double fct = 1<<sft; // (multiplication factor corresponding to sft)
//                double log_T = log2(m_maxRtt)*fct; // log2(T)*fct
//                double log_B = log2(B)*fct; // log2(B)*fct
//                double log_1e9 = log2(1e9)*fct; // log2(1e9)*fct
//                double qterm = 0;
//                double byteTerm = 0;
//                double uTerm = 0;
//                if ((qlen >> 8) > 0){
//                    int log_dt = log2apprx(dt, b, m, l); // ~log2(dt)*fct
//                    int log_qlen = log2apprx(qlen >> 8, b, m, l); // ~log2(qlen / 256)*fct
//                    qterm = pow(2, (
//                                log_dt + log_qlen + log_1e9 - log_B - 2*log_T
//                                )/fct
//                            ) * 256;
//                    // 2^((log2(dt)*fct+log2(qlen/256)*fct+log2(1e9)*fct-log2(B)*fct-2*log2(T)*fct)/fct)*256 ~= dt*qlen*1e9/(B*T^2)
//                }
//                if (m_lastPktSize[ifIndex] > 0){
//                    int byte = m_lastPktSize[ifIndex];
//                    int log_byte = log2apprx(byte, b, m, l);
//                    byteTerm = pow(2, (
//                                log_byte + log_1e9 - log_B - log_T
//                                )/fct
//                            );
//                    // 2^((log2(byte)*fct+log2(1e9)*fct-log2(B)*fct-log2(T)*fct)/fct) ~= byte*1e9 / (B*T)
//                }
//                if (m_maxRtt > dt && m_u[ifIndex] > 0){
//                    int log_T_dt = log2apprx(m_maxRtt - dt, b, m, l); // ~log2(T-dt)*fct
//                    int log_u = log2apprx(int(round(m_u[ifIndex] * 8192)), b, m, l); // ~log2(u*512)*fct
//                    uTerm = pow(2, (
//                                log_T_dt + log_u - log_T
//                                )/fct
//                            ) / 8192;
//                    // 2^((log2(T-dt)*fct+log2(u*512)*fct-log2(T)*fct)/fct)/512 = (T-dt)*u/T
//                }
//                newU = qterm+byteTerm+uTerm;
//
//                #if 0
//                /**************************
//                 * accurate calc
//                 *************************/
//                double weight_ewma = double(dt) / m_maxRtt;
//                double u;
//                if (m_lastPktSize[ifIndex] == 0)
//                    u = 0;
//                else{
//                    double txRate = m_lastPktSize[ifIndex] / double(dt); // B/ns
//                    u = (qlen / m_maxRtt + txRate) * 1e9 / B;
//                }
//                newU = m_u[ifIndex] * (1 - weight_ewma) + u * weight_ewma;
//                printf(" %lf\n", newU);
//                #endif
//
//                /************************
//                 * update PINT header
//                 ***********************/
//                uint16_t power = Pint::encode_u(newU);
//                if (power > ih->GetPower())
//                    ih->SetPower(power);
//
//                m_u[ifIndex] = newU;
//            }
//        }
//    }
//    m_txBytes[ifIndex] += p->GetSize();
//    m_lastPktSize[ifIndex] = p->GetSize();
//    m_lastPktTs[ifIndex] = Simulator::Now().GetTimeStep();
//}

int EnquserverNode::logres_shift(int b, int l){
    static int data[] = {0,0,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};
    return l - data[b];
}

int EnquserverNode::log2apprx(int x, int b, int m, int l){
    int x0 = x;
    int msb = int(log2(x)) + 1;
    if (msb > m){
        x = (x >> (msb - m) << (msb - m));
        #if 0
        x += + (1 << (msb - m - 1));
        #else
        int mask = (1 << (msb-m)) - 1;
        if ((x0 & mask) > (rand() & mask))
            x += 1<<(msb-m);
        #endif
    }
    return int(log2(x) * (1<<logres_shift(b, l)));
}

} /* namespace ns3 */
