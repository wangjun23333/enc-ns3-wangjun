#include "ns3/simulator.h"
#include "ns3/ipv4.h"
#include "ns3/packet.h"
#include "ns3/ipv4-header.h"
#include "ns3/flow-id-tag.h"
#include "ns3/boolean.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"
#include "switch-node.h"
#include "enc-net-device.h"
#include "ppp-header.h"
//#include "ns3/int-header-niux.h"
#include "../../network/utils/int-header-niux.h"
#include <cmath>

namespace ns3 {

TypeId SwitchNode::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::SwitchNode")
    .SetParent<Node> ()
    .AddConstructor<SwitchNode> ()
    .AddAttribute("EcnEnabled",
            "Enable ECN marking.",
            BooleanValue(false),
            MakeBooleanAccessor(&SwitchNode::m_ecnEnabled),
            MakeBooleanChecker())
    .AddAttribute("CcMode",
            "CC mode.",
            UintegerValue(0),
            MakeUintegerAccessor(&SwitchNode::m_ccMode),
            MakeUintegerChecker<uint32_t>())
    .AddAttribute("AckHighPrio",
            "Set high priority for ACK/NACK or not",
            UintegerValue(0),
            MakeUintegerAccessor(&SwitchNode::m_ackHighPrio),
            MakeUintegerChecker<uint32_t>())
    .AddAttribute("MaxRtt",
            "Max Rtt of the network",
            UintegerValue(9000),
            MakeUintegerAccessor(&SwitchNode::m_maxRtt),
            MakeUintegerChecker<uint32_t>())
  ;
  return tid;
}

SwitchNode::SwitchNode(uint8_t _id){
    //m_ecmpSeed = m_id;
    id = _id;

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
    for (uint32_t i = 0; i < pCnt; i++)
        max_rate[i] = 0;
}

void SwitchNode::SetMaxRate(uint8_t _port, uint64_t _max_rate) {
    max_rate[_port] = _max_rate;
    Ptr<EncNetDevice> device = DynamicCast<EncNetDevice>(m_devices[_port]);
    device->SetDataRate(_max_rate);
}

int SwitchNode::GetOutDev(Ptr<const Packet> p, CustomHeader &ch){
    // look up entries
    auto entry = m_rtTable.find(ch.dip);

    // no matching entry
    if (entry == m_rtTable.end())
        return -1;

    // entry found
    int nexthops = entry->second;

    // pick one next hop based on hash
    union {
        uint8_t u8[4+4+2+2];
        uint32_t u32[3];
    } buf;
    buf.u32[0] = ch.sip;
    buf.u32[1] = ch.dip;
    if (ch.l3Prot == 0x6)
        buf.u32[2] = ch.tcp.sport | ((uint32_t)ch.tcp.dport << 16);
    else if (ch.l3Prot == 0x11)
        buf.u32[2] = ch.udp.sport | ((uint32_t)ch.udp.dport << 16);
    else if (ch.l3Prot == 0xFC || ch.l3Prot == 0xFD)
        buf.u32[2] = ch.ack.sport | ((uint32_t)ch.ack.dport << 16);

    return nexthops;
}

void SwitchNode::CheckAndSendPfc(uint32_t inDev, uint32_t qIndex){
    Ptr<EncNetDevice> device = DynamicCast<EncNetDevice>(m_devices[inDev]);
    if (m_mmu->CheckShouldPause(inDev, qIndex)){
        m_mmu->SetPause(inDev, qIndex);
    }
}
void SwitchNode::CheckAndSendResume(uint32_t inDev, uint32_t qIndex){
    Ptr<EncNetDevice> device = DynamicCast<EncNetDevice>(m_devices[inDev]);
    if (m_mmu->CheckShouldResume(inDev, qIndex)){
        m_mmu->SetResume(inDev, qIndex);
    }
}

void SwitchNode::SendToDev(Ptr<Packet>p, CustomHeader &ch){
    int idx = GetOutDev(p, ch);
    if (idx >= 0){
        NS_ASSERT_MSG(m_devices[idx]->IsLinkUp(), "The routing table look up should return link that is up");

        // determine the qIndex
        uint32_t qIndex;
        /*if (ch.l3Prot == 0xFF || ch.l3Prot == 0xFE || (m_ackHighPrio && (ch.l3Prot == 0xFD || ch.l3Prot == 0xFC))){  //QCN or PFC or NACK, go highest priority
            qIndex = 0;
        }else{
            qIndex = (ch.l3Prot == 0x06 ? 1 : 1); // if TCP, put to queue 1
        }*/
        qIndex = 1;

        // admission control
        FlowIdTag t;
        p->PeekPacketTag(t);
        uint32_t inDev = t.GetFlowId();
        if (qIndex != 0){ //not highest priority
            if (m_mmu->CheckIngressAdmission(inDev, qIndex, p->GetSize()) && m_mmu->CheckEgressAdmission(idx, qIndex, p->GetSize())){            // Admission control
                m_mmu->UpdateIngressAdmission(inDev, qIndex, p->GetSize());
                m_mmu->UpdateEgressAdmission(idx, qIndex, p->GetSize());
            }else{
                return; // Drop
            }
            CheckAndSendPfc(inDev, qIndex);
        }
        m_bytes[inDev][idx][qIndex] += p->GetSize();
        m_devices[idx]->SwitchSend(qIndex, p, ch);
    }else
        return; // Drop
}

void SwitchNode::AddTableEntry(Ipv4Address &dstAddr, uint32_t intf_idx){
    uint32_t dip = dstAddr.Get();
    m_rtTable[dip] = intf_idx;
}

void SwitchNode::ClearTable(){
    m_rtTable.clear();
}

// This function can only be called in switch mode
bool SwitchNode::SwitchReceiveFromDevice(Ptr<NetDevice> device, Ptr<Packet> packet, CustomHeader &ch){
    SendToDev(packet, ch);
    return true;
}

void SwitchNode::SwitchNotifyDequeue(uint32_t ifIndex, uint32_t qIndex, Ptr<Packet> p){
    FlowIdTag t;
    p->PeekPacketTag(t);
    if (qIndex != 0){
        uint32_t inDev = t.GetFlowId();
        m_mmu->RemoveFromIngressAdmission(inDev, qIndex, p->GetSize());
        m_mmu->RemoveFromEgressAdmission(ifIndex, qIndex, p->GetSize());
        m_bytes[inDev][ifIndex][qIndex] -= p->GetSize();
        /*if (m_ecnEnabled){
            bool egressCongested = m_mmu->ShouldSendCN(ifIndex, qIndex);
            if (egressCongested){
                PppHeader ppp;
                Ipv4Header h;
                p->RemoveHeader(ppp);
                p->RemoveHeader(h);
                h.SetEcn((Ipv4Header::EcnType)0x03);
                p->AddHeader(h);
                p->AddHeader(ppp);
            }
        }*/
        //CheckAndSendPfc(inDev, qIndex);
        CheckAndSendResume(inDev, qIndex);
    }

    m_txBytes[ifIndex] += p->GetSize();
    m_lastPktSize[ifIndex] = p->GetSize();
    m_lastPktTs[ifIndex] = Simulator::Now().GetTimeStep();
    
    uint8_t* buf = p->GetBuffer();
    if (buf[PppHeader::GetStaticSize() + 9] == 0x06) {
        MyIntHeader *ih = (MyIntHeader*)&buf[PppHeader::GetStaticSize() + 20 + 20];
        Ptr<EncNetDevice> dev = DynamicCast<EncNetDevice>(m_devices[ifIndex]);

        uint32_t ts = m_lastPktTs[ifIndex];
        int push_rst;
        uint64_t _max_rate = max_rate[ifIndex]>>23;
        uint16_t depth = dev->GetQueue()->GetNBytesTotal();
        push_rst = ih->PushDepth(id, ifIndex, depth, ts, _max_rate);

        if (push_rst < 0) {
            uint64_t _ratio = (dev->GetDataRate().GetBitRate()*10000)/max_rate[ifIndex];
            push_rst = ih->PushRatio(id, ifIndex, _ratio, ts, _max_rate);
        }

        if (push_rst <= 0) {
            ih->PushRoute(id, ifIndex);
        }
    }
}

} /* namespace ns3 */
