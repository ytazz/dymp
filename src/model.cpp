#include <model.h>
#include <world.h>
#include <solver.h>
#include <variable.h>
#include <constraint.h>

using namespace std;

namespace dymp{;

//-------------------------------------------------------------------------------------------------

Tick::Tick(World* w, int k, real_t t){
    idx   = k;
    time  = t;
	w->ticks.push_back( unique_ptr<Tick>(this) );
}

//-------------------------------------------------------------------------------------------------

void Trajectory::Update(){
	for(auto&& key : *this){
		if(key->prev) key->hprev = key->tick->time - key->prev->tick->time;
		if(key->next) key->hnext = key->next->tick->time - key->tick->time;
	}
}

void Trajectory::Set(World* w, Model* model){
	clear();

	for(auto&& t : w->ticks)
		push_back(unique_ptr<Keypoint>(model->CreateKeypoint()));
		
	// name keypoints [node name]_i
	stringstream ss;
	for(int i = 0; i < (int)w->ticks.size(); i++){
		ss.str("");
		ss << model->name << '_' << i;
		at(i)->name = ss.str();
	}

	// then link them together
	for(int i = 0; i < (int)size(); i++){
		auto&& key = at(i);
		key->model = model;
		key->tick  = w->ticks[i].get();
		key->prev  = (i > 0		   ? (*this)[i-1].get() : NULL);
		key->next  = (i < size()-1 ? (*this)[i+1].get() : NULL);
	}

	Update();
}

void Trajectory::AddVar(Solver* solver){
	for(auto&& key : *this)
		key->AddVar(solver);
}
void Trajectory::AddCon(Solver* solver){
	for(auto&& key : *this)
		key->AddCon(solver);
}
void Trajectory::Init(){
	for(auto&& key : *this)
		key->Init();
}
void Trajectory::Prepare(){
	Update();
	for(auto&& key : *this)
		key->Prepare();
}
void Trajectory::PrepareStep(){
	for(auto&& key : *this)
		key->PrepareStep();
}
void Trajectory::Finish(){
	for(auto&& key : *this)
		key->Finish();
}

//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------

Model::Model(World* w, const string& n){
	world = w;
    w->models.push_back(unique_ptr<Model>(this));

	type  = -1;
	name  = n;

	world->ready = false;
}

Model::~Model(){
	world->ready = false;
}

void Model::Init(){

}

void Model::AddVar(){
	traj.AddVar(world->solver.get());
}
void Model::AddCon(){
	traj.AddCon(world->solver.get());
}
void Model::Prepare(){
	traj.Prepare();
}
void Model::PrepareStep(){
	traj.PrepareStep();
}
void Model::Finish(){
	traj.Finish();
}
void Model::AddKeypoints(){
	traj.Set(world, this);
}

}
