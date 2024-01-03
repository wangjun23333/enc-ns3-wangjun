#include "int-header-niux.h"

namespace ns3 {

uint32_t MyIntHeader::GetStaticSize() {
	return sizeof(hinfo)+sizeof(idInfo)+sizeof(dinfo)+sizeof(rinfo);
}

void MyIntHeader::PushRoute(uint8_t _id, uint8_t _port) {
	if (hinfo.nodeNum < idNum)
		if (rand()%4 == 0)
			iinfo[hinfo.nodeNum++].Set(_id, _port);
}

int MyIntHeader::PushDepth(uint8_t _id, uint8_t _port, uint16_t _depth, uint32_t _ts, uint8_t _maxRate) {
	_depth = _depth / qlenUnit;
	_depth = (_depth < 0xffff)? _depth : 0xffff;
	
	if (_depth <= 0) {
		return -1;
	}
	else if (hinfo.depthNum < maxNum) {
		dinfo[hinfo.depthNum++].Set(_id, _port, _depth, _ts, _maxRate);
		return 1;
	}
	else {
		uint16_t min_depth = dinfo[0].depth;
		int min_idx = 0;
		for (int i = 1; i < maxNum; ++i) {
			if (dinfo[i].depth < min_depth) {
				min_depth = dinfo[i].depth;
				min_idx = i;
			}
		}
		if (_depth > min_depth) {
			dinfo[min_idx].Set(_id, _port, _depth, _ts, _maxRate);
			return 1;
		}
		return 0;
	}
}

int MyIntHeader::PushRatio(uint8_t _id, uint8_t _port, uint16_t _ratio, uint32_t _ts, uint8_t _maxRate) {
	if (hinfo.ratioNum < maxNum) {
		rinfo[hinfo.ratioNum++].Set(_id, _port, _ratio, _ts, _maxRate);
		return 1;
	}
	else {
		uint16_t min_ratio = rinfo[0].ratio;
		int min_idx = 0;
		for (int i = 1; i < maxNum; ++i) {
			if (rinfo[i].ratio < min_ratio) {
				min_ratio = rinfo[i].ratio;
				min_idx = i;
			}
		}
		if (_ratio > min_ratio) {
			rinfo[min_idx].Set(_id, _port, _ratio, _ts, _maxRate);
			return 1;
		}
		return 0;
	}
}

void MyIntHeader::Serialize (Buffer::Iterator start) const{
	Buffer::Iterator i = start;
	i.WriteU16(hinfo.buf);
	for (int j = 0; j < idNum; ++j)
		i.WriteU16(iinfo[j].buf);
	for (int j = 0; j < maxNum; ++j)
		i.WriteU64(dinfo[j].buf);
	for (int j = 0; j < maxNum; ++j)
		i.WriteU64(rinfo[j].buf);
}

uint32_t MyIntHeader::Deserialize (Buffer::Iterator start){
	Buffer::Iterator i = start;
	hinfo.buf = i.ReadU16();
	for (int j = 0; j < idNum; ++j)
		iinfo[j].buf = i.ReadU16();
	for (int j = 0; j < maxNum; ++j)
		dinfo[j].buf = i.ReadU64();
	for (int j = 0; j < maxNum; ++j)
		rinfo[j].buf = i.ReadU64();
}

}
