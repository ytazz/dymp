#pragma once

#include <types.h>
#include <canvas.h>

namespace dymp{;

class World;
class Solver;
class Model;

/**
	tick : discrete time instants
 */
class Tick{
public:
	real_t	time;		///< time instant
	int	    idx;		///< index in tick sequence

public:
	Tick(World* w, int k, real_t t);
};

class Ticks : public std::vector< std::unique_ptr<Tick> >{
	void	AssignIndices();
public:
	void	Add   (Tick* tick);
	void	Remove(Tick* tick);
};

/**
	base class of trajectory keypoint
 */
class Keypoint{
public:
	std::string  name;
	Model*		 model;			///< reference to owner node
	Tick*		 tick;			///< time instant at which this keypoint is put
	real_t		 hprev, hnext;	///< time interval between adjacent keypoints
	Keypoint	 *prev, *next;	///< reference to adjacent keypoints	
	
public:
	/// add variables
	virtual void AddVar(Solver* solver){}
	/// add constraints
	virtual void AddCon(Solver* solver){}
	/// initialization
	virtual void Init(){}
	/// pre-processing
	virtual void Prepare(){}
	/// pre-processing
	virtual void PrepareStep(){}
	/// post-processing
	virtual void Finish(){}

	/// draw
	virtual void Draw(render::Canvas* canvas, render::Config* conf){}

	Keypoint(){
		tick  = 0;
		hprev = hnext = 0.0;
		prev  = next  = 0;
	}
};

typedef std::pair<Keypoint*, Keypoint*> KeyPair;

/**
	trajectory: a sequnce of keypoints
 */
class Trajectory : public std::vector< std::unique_ptr<Keypoint> >{
public:
	/// update length of time between each pair of keypoints
	void Update();

	/// create keypoints
	void Set(World* w, Model* model);

	/// get keypoint
	Keypoint* GetKeypoint(int idx   ){ return at(idx).get();       }
	Keypoint* GetKeypoint(Tick* tick){ return at(tick->idx).get(); }
	
	/// get segment (pair of keypoints) which includes specified time instant
	KeyPair GetSegment(real_t time){
		int idx;
		for(idx = -1; idx < (int)size()-1; idx++){
			if(time < (*this)[idx+1]->tick->time)
				break;
		}
		if(idx == -1)
			return std::make_pair(GetKeypoint(0), GetKeypoint(0));
		if(idx == size()-1)
			return std::make_pair(GetKeypoint(idx), GetKeypoint(idx));
		return std::make_pair(GetKeypoint(idx), GetKeypoint(idx+1));
	}

	void AddVar(Solver* solver);
	void AddCon(Solver* solver);
	void Init       ();
	void Prepare    ();
	void PrepareStep();
	void Finish     ();
	void Draw(render::Canvas* canvas, render::Config* conf);
};

/**
	model
 **/
class Model{
public:
    World*      world;
	std::string name;
	int         type;

	Trajectory	traj;

public:
    /// draw
	virtual void Draw(render::Canvas* canvas, render::Config* conf){}
    
    virtual void AddVar     ();
	virtual void AddCon     ();
	virtual void Init       ();
    virtual void Prepare    ();
	virtual void PrepareStep();
	virtual void Finish     ();
	
	/// create keypoint
	virtual Keypoint* CreateKeypoint() = 0;
	
	/// take snapshot at given time
	virtual void CreateSnapshot(real_t t){}

	/// draw snapshot
	virtual void DrawSnapshot(render::Canvas* canvas, render::Config* conf){}

 	void AddKeypoints();

	Model(World* w, const std::string& n);
	virtual ~Model();
};

}
