#pragma once

#include <types.h>
#include <timer.h>
#include <canvas.h>

namespace dimp3{;

class Solver;
class Tick;
class Model;

/**
	
 */
class World{
public:
	bool			ready;				///< ready flag

    std::vector< std::unique_ptr<Tick > > ticks;
	std::vector< std::unique_ptr<Model> > models;

	std::unique_ptr<Solver>	   solver;	///< internal solver
	std::unique_ptr<render::Config>  conf;	///< default draw configuration

	Timer timer;	

	int TPrepare;
	int TPrepareStep;
	int TStep;
	int TFinish;

	void Prepare    ();
	void PrepareStep();
	void Finish     ();

public:
	/// take snapshot
	void    CreateSnapshot(real_t t);

	/// does initialization
	virtual void Init();

	/// clears all objects
	virtual void Clear();

	/// resets all variables. objects are not deleted
	virtual void Reset();

	/** 
		computes one step of planning
	 **/
	virtual void Step();

	/// visualize the plan
	virtual void Draw(render::Canvas* canvas, render::Config* conf = 0);

	/// draw snapshot
	virtual void DrawSnapshot(real_t time, render::Canvas* canvas, render::Config* conf = 0);


	World();
	virtual ~World(){}
};

}
