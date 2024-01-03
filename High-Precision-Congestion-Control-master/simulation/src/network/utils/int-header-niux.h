#ifndef INT_HEADER_NIUX_H
#define INT_HEADER_NIUX_H

#include "ns3/buffer.h"
#include <stdint.h>
#include <cstdio>

namespace ns3 {

class headerInfo {
public:
	union {
		struct {
			uint16_t totalLength: 4,	// total length of this header in packet
							 nodeNum: 4,	// nodes count for routing
							 depthNum: 4,	// nodes count for depth info
							 ratioNum: 4;	// nodes count for ratio info
		};
		uint16_t buf;
	};

	headerInfo() {
		totalLength = 40;
		nodeNum = 0;
		depthNum = 0;
		ratioNum = 0;
	}
};

class idInfo {
public:
	union {
		struct {
			uint16_t id: 8,
							 port: 8;
		};
		uint16_t buf;
	};

	void Set(uint8_t _id, uint8_t _port) {
		id = _id;
		port = _port;
	}
};


class depthInfo {
public:
	union {
		struct {
			idInfo iinfo;
			uint16_t depth;
			uint32_t ts: 24,
							 maxRate: 8;
		};
		uint64_t buf;
	};
	void Set(uint8_t _id, uint8_t _port, uint16_t _depth, uint32_t _ts, uint8_t _maxRate) {
		iinfo.Set(_id, _port);
		depth = _depth;
		ts = _ts;
		maxRate = _maxRate;
	}
};

class ratioInfo {
public:
	union {
		struct {
			idInfo iinfo;
			uint16_t ratio;
			uint32_t ts: 24,
							 maxRate: 8;
		};
		uint64_t buf;
	};
	void Set(uint8_t _id, uint8_t _port, uint16_t _ratio, uint32_t _ts, uint8_t _maxRate) {
		iinfo.Set(_id, _port);
		ratio = _ratio;
		ts = _ts;
		maxRate = _maxRate;
	}
};

class MyIntHeader {
public:
	static const uint32_t idNum = 1;
	static const uint32_t maxNum = 2;

	// headerInfo: 2 Bytes
	headerInfo hinfo;
	// idInfo: 2*3 = 6 Bytes
	idInfo iinfo[idNum];
	// depthInfo: Max 8*2 = 16 Bytes
	depthInfo dinfo[maxNum];
	// ratioInfo: Max 8*2 = 16 Bytes
	ratioInfo rinfo[maxNum];

	static const uint32_t qlenUnit = 80;

	MyIntHeader();
	static uint32_t GetStaticSize();
	void PushRoute(uint8_t _id, uint8_t _port);
	int PushDepth(uint8_t _id, uint8_t _port, uint16_t _depth, uint32_t _ts, uint8_t _maxRate);
	int PushRatio(uint8_t _id, uint8_t _port, uint16_t _ratio, uint32_t _ts, uint8_t _maxRate);
	void Serialize (Buffer::Iterator start) const;
	uint32_t Deserialize (Buffer::Iterator start);
};

}

#endif /* INT_HEADER_H */
