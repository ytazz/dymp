#include <dymp/world.h>
#include <dymp/solver.h>
#include <dymp/variable.h>
#include <dymp/constraint.h>
#include <dymp/model.h>

using namespace std;

namespace dymp{;

World::World(){
	solver = unique_ptr<Solver>(new Solver());
    conf   = unique_ptr<render::Config>(new render::Config());
}

void World::Init(){
	//omp_set_num_threads(1);

	//solver.numthread = param.numthread;
	solver->Clear();

	// create keypoints of trajectory nodes
    for(auto&& m : models)
        m->AddKeypoints();

	// register variables to solver
	for(auto&& m : models)
        m->AddVar();

	// register constraints to solver
	for(auto&& m : models)
        m->AddCon();

	// do extra initialization
	for(auto&& m : models)
        m->Init();

	solver->Init();

	ready = true;
}

void World::Clear(){
	models.clear();
    solver->Clear();
	ready = false;
}

void World::Reset(){
	solver->Reset();
}

void World::Prepare(){
	for(auto&& m : models)
        m->Prepare();
}

void World::PrepareStep(){
	for(auto&& m : models)
        m->PrepareStep();
}

void World::Finish(){
	for(auto&& m : models)
        m->Finish();
}

void World::Step(){
	if(!ready)
		Init();

	timer.CountUS();
	Prepare();
	TPrepare = timer.CountUS();
	
	timer.CountUS();
	PrepareStep();
	TPrepareStep = timer.CountUS();
	
	timer.CountUS();
	solver->Step();
	TStep = timer.CountUS();

	timer.CountUS();
	Finish();
	TFinish = timer.CountUS();

	if(solver->param.verbose){
		//DSTR << " tpre1: " << TPrepare;
		//DSTR << " tpre2: " << TPrepareStep;
		//DSTR << " tstp: "  << TStep;
		//DSTR << " tfin: "  << TFinish << endl;
	}

}

void World::Draw(render::Canvas* canvas, render::Config* _conf){
	if(!_conf)
		_conf = conf.get();

	if(conf->Set(canvas, render::Item::GlobalAxis, 0)){
		float l = conf->Scale(render::Item::GlobalAxis, 0);
		canvas->Line(vec3_t::Zero(), l*vec3_t::UnitX());
		canvas->Line(vec3_t::Zero(), l*vec3_t::UnitY());
		canvas->Line(vec3_t::Zero(), l*vec3_t::UnitZ());
	}

	for(auto&& m : models)
        m->Draw(canvas, _conf);
	
	for(auto&& tick : ticks)
		DrawSnapshot(tick->time, canvas, _conf);
}

void World::CreateSnapshot(real_t t){
	for(auto&& m : models)
        m->CreateSnapshot(t);
}

void World::DrawSnapshot(real_t time, render::Canvas* canvas, render::Config* _conf){
	if(!_conf)
		_conf = conf.get();

	for(auto&& m : models)
        m->CreateSnapshot(time);
    for(auto&& m : models)
        m->DrawSnapshot(canvas, _conf);
}

}
