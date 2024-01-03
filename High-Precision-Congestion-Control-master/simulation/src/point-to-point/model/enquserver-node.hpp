//
//  enquserver-node.hpp
//  NS3-cc
//
//  Created by 王钧 on 2023/9/26.
//

#ifndef ENQUSERVER_NODE_H
#define ENQUSERVER_NODE_H

#include <unordered_map>
#include <ns3/node.h>
#include "qbb-net-device.h"
#include "switch-mmu.h"
#include "pint.h"
#include <vector>

namespace ns3 {

class Packet;

class EnquserverNode : public Node{
    static const uint32_t pCnt = 257;    // Number of ports used
    static const uint32_t qCnt = 8;    // Number of queues/priorities used
    static const uint32_t eqsPort = 257;    // 出入口服务器发送携带信息数据包的源端口号
    static const uint32_t eqsIp = 010010;    // 出入口服务器发送携带信息数据包的源IP
    uint32_t m_ecmpSeed;
    
//    std::unordered_map<uint32_t, std::vector<std::pair<int, int>>> m_sharedTable; 创建一个名为m_sharedTable的数据结构，用于存储共享链路表信息，其中使用32位整数作为键，用来表示广域网节点ID。与广域网节点ID相关联的值为一个向量，该向量存储每个发送方主机的IP和端口号。
//     monitor of PFC
//    uint32_t m_bytes[pCnt][pCnt][qCnt]; // m_bytes[inDev][outDev][qidx] is the bytes from inDev enqueued for outDev at qidx
//
//    uint64_t m_txBytes[pCnt]; // counter of tx bytes
//
//    uint32_t m_lastPktSize[pCnt];
//    uint64_t m_lastPktTs[pCnt]; // ns
//    double m_u[pCnt];
    
    std::unordered_map<int32_t, int> m_routerMap;
    
    struct flowInfo{
        uint32_t sip;
        uint32_t dip;
        uint16_t sport;
        uint16_t dport;
        // 实现比较操作符
        bool operator==(const flowInfo& other) const {
            return sip == other.sip && dip == other.dip &&
                   sport == other.sport && dport == other.dport;
        }

        // 实现小于操作符，用于 set 的排序
        bool operator<(const flowInfo& other) const {
            return std::tie(sip, dip, sport, dport) < std::tie(other.sip, other.dip, other.sport, other.dport);
        }
    };
    
    
    struct m_sharedTableEntry{
        uint16_t rid;
        uint16_t port;
        std::vector<flowInfo> flowInfos;
    };
    
    struct relatedSenderHeaderInfo{
        flowInfo fInfo;
        std::vector<std::pair<uint16_t, uint16_t>> rIdAndPort;
    };
    
    std::vector<m_sharedTableEntry> m_sharedTable;
//
//
//    struct HeaderLinkInfo{
//        uint32_t eqsPort = 257;
//        uint32_t eqsIp = 010010;
//        uint32_t dip;
//        uint32_t dport;
//        std::pair<int, int> rate_ratio1; // 使用 std::pair
//        std::pair<int, int> rate_ratio2;
//        std::pair<int, int> rate_ratio3;
//        std::pair<int, int> rate_ratio4;
//        std::tuple<int, int, int> qLen1;
//        std::tuple<int, int, int> qLen2;
//        std::tuple<int, int, int> qLen3;
//        std::tuple<int, int, int> qLen4;
//    };
protected:
    bool m_ecnEnabled;
    uint32_t m_ccMode;
    uint64_t m_maxRtt;

    uint32_t m_ackHighPrio; // set high priority for ACK/NACK

private:
    int GetOutDev(Ptr<const Packet>p, CustomHeader &ch);
    void SendToDev(Ptr<Packet>p, CustomHeader &ch);
    static uint32_t EcmpHash(const uint8_t* key, size_t len, uint32_t seed);
    void CheckAndSendPfc(uint32_t inDev, uint32_t qIndex);
    void CheckAndSendResume(uint32_t inDev, uint32_t qIndex);
    void GetShareTable(Ptr<const Packet>p, CustomHeader &ch);//获取共享链路表的函数，参数为数据包包头的广域网节点ID，目的地址端口号
//    void MatchSharedTableSendToRelatedSender(Ptr<Packet>p, CustomHeader &ch);
    //对携带链路信息的数据包中的信息和共享链路表进行查找匹配，返回HeaderLinkInfo结构体类型中的数据
    
public:
    Ptr<SwitchMmu> m_mmu;

    static TypeId GetTypeId (void);
    EnquserverNode();
    void SetEcmpSeed(uint32_t seed);
    void AddTableEntry(Ipv4Address &dstAddr, uint32_t intf_idx);
    void ClearTable();
//    bool SwitchReceiveFromDevice(Ptr<NetDevice> device, Ptr<Packet> packet, CustomHeader &ch);
    void MatchSharedTableSendToRelatedSender(Ptr<NetDevice> device, Ptr<Packet>p, CustomHeader &ch);

    // for approximate calc in PINT
    int logres_shift(int b, int l);
    int log2apprx(int x, int b, int m, int l); // given x of at most b bits, use most significant m bits of x, calc the result in l bits
};

} /* namespace ns3 */

#endif /* ENQUSERVER_NODE_H */
