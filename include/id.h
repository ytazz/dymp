#pragma once

#include <string>

namespace dymp{;

/// variable identifiers
struct VarTag{
	enum{
		Any = 0,
		CentroidPosT,
		CentroidPosR,
		CentroidVelT,
		CentroidVelR,
		CentroidMomentum,
		CentroidTime,
		CentroidDuration,
		CentroidEndPos,
		CentroidEndVel,
		CentroidEndAcc,
		CentroidEndStiff,
		CentroidEndCmp,
		CentroidEndForce,
		CentroidEndMoment,
		WholebodyPosT,
		WholebodyPosR,
		WholebodyVelT,
		WholebodyVelR,
		WholebodyAccT,
		WholebodyAccR,
        WholebodyMomentum,
		WholebodyJointPos,
		WholebodyJointVel,
		WholebodyJointAcc,
		WholebodyJointJerk,
		WholebodyForceT,
		WholebodyForceR,
		NumTypes,
	};
};

/// constraint identifiers
struct ConTag{
	enum{
		Any = 0,
		CentroidPosT,
		CentroidPosR,
		CentroidVelT,
		CentroidVelR,
		CentroidMomentum,
		CentroidTime,
        CentroidEndPosT,
		CentroidEndPosR,
		CentroidDesPosT,
		CentroidDesPosR,
		CentroidDesVelT,
		CentroidDesVelR,
		CentroidDesMomentum,
		CentroidDesTime,
        CentroidDesDuration,
        CentroidDesEndPosT,
		CentroidDesEndVelT,
		CentroidDesEndPosR,
		CentroidDesEndVelR,
		CentroidDesEndStiff,
		CentroidDesEndCmp,
		CentroidDesEndForce,
		CentroidDesEndMoment,
		CentroidDurationRange,
		CentroidEndPosRange,
		CentroidEndStiffRange,
		CentroidEndContact,
		CentroidEndFriction,
		CentroidEndMomentRange,
		WholebodyPosT,
		WholebodyPosR,
		WholebodyVelT,
		WholebodyVelR,
		WholebodyAccT,
		WholebodyAccR,
		WholebodyJointPos,
		WholebodyJointVel,
		WholebodyJointAcc,
		WholebodyJointJerk,
		WholebodyForceT,
		WholebodyForceR,
		WholebodyLimit,
		WholebodyContactPosT,
		WholebodyContactPosR,
		WholebodyContactVelT,
		WholebodyContactVelR,
		WholebodyNormalForce,
		WholebodyFrictionForce,
		WholebodyMoment,
		WholebodyMomentum,
		NumTypes,
	};
};

extern const char* VarNames[VarTag::NumTypes];
extern const char* ConNames[ConTag::NumTypes];

class ID{
public:
	int          tag;
	void*        owner;
	std::string  name;

	int Match(ID* id){
		if(tag == 0)
			return 1;
		if(tag != id->tag)
			return 0;
		if(!owner)
			return 2;
		if(owner != id->owner)
			return 0;
		return 3;
	}

	ID(int _tag = 0, void* _owner = 0, std::string _name = ""){
		tag    = _tag;
		owner  = _owner;
		name   = _name ;
	}
};

}
