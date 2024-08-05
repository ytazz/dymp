#include <app.h>
#include <centroid.h>
#include <solver.h>
#include <util.h>
#include <rollpitchyaw.h>

/**
 centroidal trajectory planning with unscheduled contact
*/

namespace dymp{;

class MyApp : public App{
public:
    struct Scene{
        enum{
            Flat,
            Gap,
            GapWithRail,
            GapWithWall,
            Steps,
            Stairs,
        };
    };
    struct Task{
        enum{
            Travel,
            LongJump,
            Backflip,
        };
    };
    struct Gait{
        enum{
            WalkWithDoubleSupport,
            WalkWithoutDoubleSupport,
            Run,
            TrotWithQuadSupport,
            TrotWithoutQuadSupport,
            Pace,
            Gallop,
        };
    };
    struct Robot{
        enum{
            Biped,
            Humanoid,
            Quadruped,
        };
    };
    struct EndConfig{
        vec3_t basePos  ;
        vec3_t posOrigin;
        vec3_t posMin   ;
        vec3_t posMax   ;
        real_t stiffMax ;
        vec2_t cmpOffset;
    };

    Centroid*	       centroid;
    FILE*              fileDuration;
    FILE*              fileCost;
    vector<EndConfig>  endConf;
        
    real_t comHeight;
    vec3_t startPos;
    vec3_t startOri;
    vec3_t goalPos;
    vec3_t goalOri;
    int    N;
    real_t dt;
    real_t goalTime;
    int    robotSelect;
    int    sceneSelect;
    int    taskSelect;
    int    gaitSelect;

public:
    virtual void BuildScene(){
        robotSelect = Robot::Biped;
        //robotSelect = Robot::Humanoid;
        //robotSelect = Robot::Quadruped;
        sceneSelect = Scene::Flat;
        //sceneSelect = Scene::Gap;
        //sceneSelect = Scene::Stairs;
        //sceneSelect = Scene::Steps;
        taskSelect = Task::Travel;
        //taskSelect = Task::LongJump;
        //taskSelect = Task::Backflip;
        //taskSelect = Task::Turn;
        //gaitSelect = Gait::WalkWithDoubleSupport;
        gaitSelect = Gait::Run;
        //gaitSelect = Gait::TrotWithQuadSupport;
        //gaitSelect = Gait::TrotWithoutQuadSupport;
        //gaitSelect = Gait::Pace;
        //gaitSelect = Gait::Gallop;

        if(robotSelect == Robot::Biped){
            comHeight = 0.7;
            endConf.resize(2);
            endConf[0].basePos   = vec3_t( 0.0, -0.15/2.0,  0.0);
            endConf[1].basePos   = vec3_t( 0.0,  0.15/2.0,  0.0);
            endConf[0].posOrigin = vec3_t( 0.0,  0.0, -comHeight);
            endConf[1].posOrigin = vec3_t( 0.0,  0.0, -comHeight);
            endConf[0].posMin    = vec3_t(-0.3, -0.0, -comHeight-0.1);
            endConf[1].posMin    = vec3_t(-0.3, -0.0, -comHeight-0.1);
            endConf[0].posMax    = vec3_t( 0.3,  0.0, -comHeight+0.3);
            endConf[1].posMax    = vec3_t( 0.3,  0.0, -comHeight+0.3);
            endConf[0].stiffMax  = 50.0;
            endConf[1].stiffMax  = 50.0;
            endConf[0].cmpOffset = -0.0*vec2_t(endConf[0].basePos.x(), endConf[0].basePos.y());
            endConf[1].cmpOffset = -0.0*vec2_t(endConf[1].basePos.x(), endConf[1].basePos.y());
        }
        if(robotSelect == Robot::Humanoid){
            comHeight = 0.7;
            endConf.resize(4);
            endConf[0].basePos   = vec3_t( 0.0, -0.20/2.0,  0.0);
            endConf[1].basePos   = vec3_t( 0.0,  0.20/2.0,  0.0);
            endConf[2].basePos   = vec3_t( 0.0, -0.25/2.0,  0.5);
            endConf[3].basePos   = vec3_t( 0.0,  0.25/2.0,  0.5);
            endConf[0].posOrigin = vec3_t( 0.0,  0.0, -0.7);
            endConf[1].posOrigin = vec3_t( 0.0,  0.0, -0.7);
            endConf[2].posOrigin = vec3_t( 0.0, -0.1, -0.5);
            endConf[3].posOrigin = vec3_t( 0.0,  0.1, -0.5);
            endConf[0].posMin    = vec3_t(-0.3, -0.05, -0.75);
            endConf[1].posMin    = vec3_t(-0.3, -0.05, -0.75);
            endConf[2].posMin    = vec3_t(-0.2, -0.3, -0.6);
            endConf[3].posMin    = vec3_t(-0.2,  0.0, -0.6);
            endConf[0].posMax    = vec3_t( 0.3,  0.05, -0.55);
            endConf[1].posMax    = vec3_t( 0.3,  0.05, -0.55);
            endConf[2].posMax    = vec3_t( 0.2,  0.0,  0.5);
            endConf[3].posMax    = vec3_t( 0.2,  0.3,  0.5);
            endConf[0].stiffMax  = 50.0;
            endConf[1].stiffMax  = 50.0;
            endConf[2].stiffMax  = 50.0;
            endConf[3].stiffMax  = 50.0;
        }
        if(robotSelect == Robot::Quadruped){
            comHeight = 0.5;
            endConf.resize(4);
            endConf[0].basePos   = vec3_t( 0.7/2.0, -0.4/2.0,  0.0);
            endConf[1].basePos   = vec3_t( 0.7/2.0,  0.4/2.0,  0.0);
            endConf[2].basePos   = vec3_t(-0.7/2.0, -0.4/2.0,  0.0);
            endConf[3].basePos   = vec3_t(-0.7/2.0,  0.4/2.0,  0.0);
            endConf[0].posOrigin = vec3_t( 0.0,  0.0, -0.5);
            endConf[1].posOrigin = vec3_t( 0.0,  0.0, -0.5);
            endConf[2].posOrigin = vec3_t( 0.0,  0.0, -0.5);
            endConf[3].posOrigin = vec3_t( 0.0,  0.0, -0.5);
            endConf[0].posMin    = vec3_t(-0.30, -0.15, -0.6);
            endConf[1].posMin    = vec3_t(-0.30, -0.15, -0.6);
            endConf[2].posMin    = vec3_t(-0.30, -0.15, -0.6);
            endConf[3].posMin    = vec3_t(-0.30, -0.15, -0.6);
            endConf[0].posMax    = vec3_t( 0.30,  0.15, -0.4);
            endConf[1].posMax    = vec3_t( 0.30,  0.15, -0.4);
            endConf[2].posMax    = vec3_t( 0.30,  0.15, -0.4);
            endConf[3].posMax    = vec3_t( 0.30,  0.15, -0.4);
            endConf[0].stiffMax  = 50.0;
            endConf[1].stiffMax  = 50.0;
            endConf[2].stiffMax  = 50.0;
            endConf[3].stiffMax  = 50.0;
		    endConf[0].cmpOffset = vec2_t(-1.0*endConf[0].basePos.x(), -1.0*endConf[0].basePos.y());
            endConf[1].cmpOffset = vec2_t(-1.0*endConf[1].basePos.x(), -1.0*endConf[1].basePos.y());
            endConf[2].cmpOffset = vec2_t(-1.0*endConf[2].basePos.x(), -1.0*endConf[2].basePos.y());
            endConf[3].cmpOffset = vec2_t(-1.0*endConf[3].basePos.x(), -1.0*endConf[3].basePos.y());
        }
        centroid = new Centroid(world.get(), "centroid");
        
	    centroid->param.g = 9.8;
	    centroid->param.m = 44.0;
        centroid->param.I(0,0) = 5.0;
        centroid->param.I(1,1) = 5.0;
        centroid->param.I(2,2) = 5.0;
        centroid->param.swingHeight = 0.15;
        centroid->param.mu = 1.0;

        centroid->param.durationMin = 0.3;
        centroid->param.durationMax = 0.7;

        centroid->param.contactMargin = 0.0;
        centroid->param.complWeight   = 1000.0;

        centroid->param.enableRotation   = true;
        centroid->param.rotationResolution = 10;
        //centroid->param.endInterpolation = DiMP::Centroid::EndInterpolation::Local;
        centroid->param.endInterpolation = Centroid::EndInterpolation::Global;
        centroid->param.endWrenchParametrization = Centroid::EndWrenchParametrization::Stiffness;
        //centroid->param.endWrenchParametrization = DiMP::Centroid::EndWrenchParametrization::Direct;
        
        // create geometry
        centroid->param.bodyRangeMin = vec3_t(-0.05, -0.05,  0.0);
        centroid->param.bodyRangeMax = vec3_t( 0.05,  0.05,  0.5);
		
	    Centroid::Face face;
        // flat ground
        if( sceneSelect == Scene::Flat ){
            real_t r = 0.0;
            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);
        }
        // gap
        if(sceneSelect == Scene::Gap){
            real_t r = 0.0;
            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);
        
            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);
        }
        // gap with rail
        if(sceneSelect == Scene::GapWithRail){
            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);
        
            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);
        
            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);

            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);
        }
        // gap with wall
        if(sceneSelect == Scene::GapWithWall){
            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);
        
            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);
        
            face.normal = vec3_t(0.0, 1.0, 0.0);
            centroid->faces.push_back(face);

            face.normal = vec3_t(0.0, -1.0, 0.0);
            centroid->faces.push_back(face);
        }
        if(sceneSelect == Scene::Stairs){
            real_t r = 0.01;
            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);
        
            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);

            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);

            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);

            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);
        }
	    if(sceneSelect == Scene::Steps){
            real_t r = 0.01;
            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);
        
            face.normal = vec3_t(0.0, 0.0, 1.0);
            centroid->faces.push_back(face);
        }
		
        const int infi = numeric_limits<int>::max();
        int nend = endConf.size();
		
        if(taskSelect == Task::Travel){
            goalTime = 8.0;
            
            if(sceneSelect == Scene::Flat){
                startPos = vec3_t(0.0, 0.0, comHeight);
		        startOri = vec3_t(0.0, 0.0, 0.0);
                if(robotSelect == Robot::Biped){
		            //goalPos  = vec3_t(3.0, 0.0, comHeight);
                    goalPos  = vec3_t(5.0, 0.0, comHeight);
                    goalOri  = vec3_t(0.0, 0.0, deg_to_rad(0.0));
		            //goalOri  = vec3_t(0.0, 0.0, Rad(180.0));
                }
                if(robotSelect == Robot::Humanoid){
		            goalPos  = vec3_t(3.0, 3.0, comHeight);
		            goalOri  = vec3_t(0.0, 0.0, deg_to_rad(90.0));
                }
                if(robotSelect == Robot::Quadruped){
		            goalPos  = vec3_t(3.0, 0.0, comHeight);
		            goalOri  = vec3_t(0.0, 0.0, 0.0);
                }

                if(robotSelect == Robot::Biped){
                    int ndiv_ssp = 1;
                    int ndiv_dsp = 1;
                    int ndiv_fp  = 1;
                    if(gaitSelect == Gait::WalkWithDoubleSupport){
                        centroid->phases = {
                            {{ 0, 0}, ndiv_dsp},
                            {{ 0,-1}, ndiv_ssp},
                            {{ 0, 0}, ndiv_dsp},
                            {{-1, 0}, ndiv_ssp},
                            {{ 0, 0}, ndiv_dsp},
                            {{ 0,-1}, ndiv_ssp},
                            {{ 0, 0}, ndiv_dsp},
                            {{-1, 0}, ndiv_ssp},
                            {{ 0, 0}, ndiv_dsp},
                            {{ 0,-1}, ndiv_ssp},
                            {{ 0, 0}, ndiv_dsp},
                            {{-1, 0}, ndiv_ssp},
                            {{ 0, 0}, ndiv_dsp},
                            {{ 0,-1}, ndiv_ssp},
                            {{ 0, 0}, ndiv_dsp},
                            {{-1, 0}, ndiv_ssp},
                            {{ 0, 0}, ndiv_dsp},
                            {{ 0,-1}, ndiv_ssp},
                            {{ 0, 0}, ndiv_dsp},
                            {{-1, 0}, ndiv_ssp},
                            {{ 0, 0}, ndiv_dsp},
                            {{ 0,-1}, ndiv_ssp},
                            {{ 0, 0}, ndiv_dsp},
                            {{-1, 0}, ndiv_ssp},
                            {{ 0, 0}, 1}
                        };
                    }
                    if(gaitSelect == Gait::WalkWithoutDoubleSupport){
                        centroid->phases = {
                            {{ 0,  0}, ndiv_dsp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1,  0}, ndiv_ssp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1,  0}, ndiv_ssp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1,  0}, ndiv_ssp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1,  0}, ndiv_ssp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1,  0}, ndiv_ssp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1,  0}, ndiv_ssp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1,  0}, ndiv_ssp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1,  0}, ndiv_ssp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1,  0}, ndiv_ssp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1,  0}, ndiv_ssp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1,  0}, ndiv_ssp},
                            {{ 0, -1}, ndiv_ssp},
                            {{ 0,  0}, ndiv_dsp}
                        };
                    }
                    if(gaitSelect == Gait::Run){
                        centroid->phases = {
                            {{ 0,  0}, ndiv_dsp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1, -1}, ndiv_fp},
                            {{-1,  0}, ndiv_ssp},
                            {{-1, -1}, ndiv_fp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1, -1}, ndiv_fp},
                            {{-1,  0}, ndiv_ssp},
                            {{-1, -1}, ndiv_fp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1, -1}, ndiv_fp},
                            {{-1,  0}, ndiv_ssp},
                            {{-1, -1}, ndiv_fp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1, -1}, ndiv_fp},
                            {{-1,  0}, ndiv_ssp},
                            {{-1, -1}, ndiv_fp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1, -1}, ndiv_fp},
                            {{-1,  0}, ndiv_ssp},
                            {{-1, -1}, ndiv_fp},
                            {{ 0, -1}, ndiv_ssp},
                            {{-1, -1}, ndiv_fp},
                            {{-1,  0}, ndiv_ssp},
                            {{ 0,  0}, ndiv_dsp}
                        };
                    }
                }
                if(robotSelect == Robot::Quadruped){
                    int ndiv_qsp = 1;
                    int ndiv_dsp = 3;
                    int ndiv_ssp = 1;
                    int ndiv_fp  = 1;
                    if(gaitSelect == Gait::TrotWithQuadSupport){
                        centroid->phases = {
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp}
                        };
                    }
                    if(gaitSelect == Gait::TrotWithoutQuadSupport){
                        centroid->phases = {
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{-1,  0,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp}
                        };  
                    }
                    if(gaitSelect == Gait::Pace){
                        centroid->phases = {
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1,  0, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1,  0, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1,  0, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1,  0, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1,  0, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{ 0, -1,  0, -1}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1,  0, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp}
                        };
                    }
                    if(gaitSelect == Gait::Gallop){
                        centroid->phases = {
                            {{ 0,  0,  0,  0}, ndiv_qsp},
                            {{-1, -1,  0,  0}, ndiv_dsp},
                            {{-1, -1, -1,  0}, ndiv_ssp},
                            {{-1, -1, -1, -1}, ndiv_fp },
                            {{ 0, -1, -1, -1}, ndiv_ssp},
                            {{-1,  0, -1, -1}, ndiv_ssp},
                            {{-1, -1,  0, -1}, ndiv_ssp},
                            {{-1, -1, -1,  0}, ndiv_ssp},
                            {{-1, -1, -1, -1}, ndiv_fp },
                            {{ 0, -1, -1, -1}, ndiv_ssp},
                            {{-1,  0, -1, -1}, ndiv_ssp},
                            {{-1, -1,  0, -1}, ndiv_ssp},
                            {{-1, -1, -1,  0}, ndiv_ssp},
                            {{-1, -1, -1, -1}, ndiv_fp },
                            {{ 0, -1, -1, -1}, ndiv_ssp},
                            {{-1,  0, -1, -1}, ndiv_ssp},
                            {{-1, -1,  0, -1}, ndiv_ssp},
                            {{-1, -1, -1,  0}, ndiv_ssp},
                            {{-1, -1, -1, -1}, ndiv_fp },
                            {{ 0, -1, -1, -1}, ndiv_ssp},
                            {{-1,  0, -1, -1}, ndiv_ssp},
                            {{-1, -1,  0, -1}, ndiv_ssp},
                            {{-1, -1, -1,  0}, ndiv_ssp},
                            {{ 0, -1, -1,  0}, ndiv_dsp},
                            {{ 0,  0,  0,  0}, ndiv_qsp}
                        };
                    }
                }
            }
            if( sceneSelect == Scene::Gap ||
                sceneSelect == Scene::GapWithRail ||
                sceneSelect == Scene::GapWithWall){

                startPos = vec3_t(0.0, 0.0, comHeight);
		        startOri = vec3_t(0.0, 0.0, 0.0);
                if(robotSelect == Robot::Biped){
		            goalPos  = vec3_t(3.0, 0.0, comHeight);
		            goalOri  = vec3_t(0.0, 0.0, 0.0);
                }
                if(robotSelect == Robot::Humanoid){
		            goalPos  = vec3_t(3.0, 0.0, comHeight);
		            goalOri  = vec3_t(0.0, 0.0, 0.0);
                }
                if(robotSelect == Robot::Quadruped){
		            goalPos  = vec3_t(1.5, 0.0, comHeight);
		            goalOri  = vec3_t(0.0, 0.0, 0.0);
                }
            }
            if(sceneSelect == Scene::Stairs){
                startPos = vec3_t(0.0, 0.0, comHeight);
		        startOri = vec3_t(0.0, 0.0, 0.0);
		        goalPos  = vec3_t(1.6, 0.0, comHeight + 0.4);
                goalOri  = vec3_t(0.0, 0.0, 0.0);
            }

            int N = centroid->phases.size()-1;
            real_t dt = goalTime/N;
            centroid->waypoints.resize(N+1);
            {
                Centroid::Waypoint& wp = centroid->waypoints[0];
                wp.value  = Centroid::Waypoint::Value (0.0, dt, startPos, startOri, zero3, zero3);
                wp.weight = Centroid::Waypoint::Weight(10.0, 1.0, 10.0*one3, 10.0*one3, 10.0*one3, 10.0*one3);
                wp.ends.resize(nend);
                for(int iend = 0; iend < nend; iend++){
		            wp.ends[iend].value  = Centroid::Waypoint::End::Value (startPos + endConf[iend].basePos + endConf[iend].posOrigin, zero3, zero3, zero3, infi);
                    wp.ends[iend].weight = Centroid::Waypoint::End::Weight(10.0*one3, 10.0*one3, 10.0*one3, 10.0*one3, 1.0, 1.0*one2, 1*one3);
                }
            }
            {
                Centroid::Waypoint& wp = centroid->waypoints[1];
                wp.weight.time  = 1.0;
                wp.weight.pos_t = 1*one3;
                wp.weight.pos_r = 1*one3;
                wp.weight.vel_t = 10*one3;
                wp.weight.L     = 10*one3;
                wp.ends.resize(nend);
                for(int iend = 0; iend < nend; iend++){
                    wp.ends[iend].weight.pos_t  = 1*one3;
                    wp.ends[iend].weight.pos_r  = 1*one3;
                    wp.ends[iend].weight.vel_t  = 1.0*one3;
                    wp.ends[iend].weight.vel_r  = 1.0*one3;
                }
            }
            {
                Centroid::Waypoint& wp = centroid->waypoints[N-1];
                wp.ends.resize(nend);
                wp.weight.time  = 1.0;
                wp.weight.pos_t = 1*one3;
                wp.weight.pos_r = 1*one3;
                wp.weight.vel_t = 10*one3;
                wp.weight.L     = 10*one3;
                for(int iend = 0; iend < nend; iend++){
                    wp.ends[iend].weight.pos_t  = 1*one3;
                    wp.ends[iend].weight.pos_r  = 1*one3;
                    wp.ends[iend].weight.vel_t  = 1*one3;
                    wp.ends[iend].weight.vel_r  = 1*one3;
                }
            }
            {
                Centroid::Waypoint& wp = centroid->waypoints[N];
                wp.value  = Centroid::Waypoint::Value (dt*N, dt, goalPos, goalOri, zero3, zero3);
                wp.weight = Centroid::Waypoint::Weight(1.0, 1.0, 10.0*one3, 10.0*one3, 10.0*one3, 10.0*one3);
                wp.ends.resize(nend);
                for(int iend = 0; iend < nend; iend++){
		            wp.ends[iend].value  = Centroid::Waypoint::End::Value (goalPos + FromRollPitchYaw(goalOri)*(endConf[iend].basePos + endConf[iend].posOrigin), goalOri, zero3, zero3, infi);
                    wp.ends[iend].weight = Centroid::Waypoint::End::Weight(10.0*one3, 10.0*one3, 10.0*one3, 10.0*one3, 1.0, 1.0*one2, 1.0*one3);
                }
            }
        }
        if(taskSelect == Task::LongJump){
            goalTime = 1.5;
            if(sceneSelect == Scene::Flat){
                startPos = vec3_t(0.0, 0.0, comHeight);
		        startOri = vec3_t(0.0, 0.0, 0.0);
                if(robotSelect == Robot::Biped){
		            goalPos  = vec3_t(1.0, 0.0, comHeight);
		            goalOri  = vec3_t(0.0, 0.0, deg_to_rad(0.0));
                    int ndiv_dsp = 1;
                    centroid->phases = {
                        {{ 0,  0}, ndiv_dsp},
                        {{ 0,  0}, ndiv_dsp},
                        {{-1, -1}, 1},
                        {{ 0,  0}, ndiv_dsp},
                        {{ 0,  0}, ndiv_dsp},
                        {{ 0,  0}, 1}
                    };
                }
            }
            int N = centroid->phases.size()-1;
            dt = goalTime/((real_t)N);
		    centroid->waypoints.resize(N+1);
            {
                Centroid::Waypoint& wp = centroid->waypoints[0];
                wp.value  = Centroid::Waypoint::Value (0.0, dt, startPos, startOri, zero3, zero3);
                wp.weight = Centroid::Waypoint::Weight(10.0, 1, 10.0*one3, 10.0*one3, 10.0*one3, 10.0*one3);
                wp.ends.resize(nend);
                for(int iend = 0; iend < nend; iend++){
		            wp.ends[iend].value  = Centroid::Waypoint::End::Value (startPos + endConf[iend].basePos + endConf[iend].posOrigin, zero3, zero3, zero3, infi);
                    wp.ends[iend].weight = Centroid::Waypoint::End::Weight(10.0*one3, 10.0*one3, 10.0*one3, 10.0*one3, 1.0, 1.0*one2, 1.0*one3);
                }
            }
            {
                Centroid::Waypoint& wp = centroid->waypoints[1];
                wp.weight.time  = 1;
                wp.weight.pos_t = 1*one3;
                wp.weight.pos_r = 1*one3;
                wp.weight.vel_t = 1*one3;
                wp.weight.L     = 1*one3;
                wp.ends.resize(nend);
                for(int iend = 0; iend < nend; iend++){
                    wp.ends[iend].weight.pos_t  = 1*one3;
                    wp.ends[iend].weight.vel_t  = 1*one3;
                }
            }
            {
                Centroid::Waypoint& wp = centroid->waypoints[N-1];
                wp.ends.resize(nend);
                wp.weight.time  = 1;
                wp.weight.pos_t = 1*one3;
                wp.weight.pos_r = 1*one3;
                wp.weight.vel_t = 1*one3;
                wp.weight.L     = 1*one3;
                for(int iend = 0; iend < nend; iend++){
                    wp.ends[iend].weight.pos_t  = 1*one3;
                    wp.ends[iend].weight.vel_t  = 1*one3;
                }
            }
            {
                Centroid::Waypoint& wp = centroid->waypoints[N];
                wp.value  = Centroid::Waypoint::Value (dt*N, dt, goalPos, goalOri, zero3, zero3);
                wp.weight = Centroid::Waypoint::Weight(1, 1, 10.0*one3, 10.0*one3, 10.0*one3, 10.0*one3);
                wp.ends.resize(nend);
                for(int iend = 0; iend < nend; iend++){
		            wp.ends[iend].value  = Centroid::Waypoint::End::Value (goalPos + endConf[iend].basePos + endConf[iend].posOrigin, zero3, zero3, zero3, infi);
                    wp.ends[iend].weight = Centroid::Waypoint::End::Weight(10.0*one3, 10.0*one3, 10.0*one3, 10.0*one3, 1.0, 1.0*one2, 1.0*one3);
                }
            }
        }
        if(taskSelect == Task::Backflip){
		    if(sceneSelect == Scene::Steps){
                startPos = vec3_t( 0.25, 0.0, comHeight + 0.3);
                startOri = vec3_t( 0.0, 0.0, 0.0);
                goalPos  = vec3_t(-0.25, 0.0, comHeight);
                goalOri  = vec3_t(0.0, deg_to_rad(-360.0), 0.0);
                if(robotSelect == Robot::Biped){
                    int ndiv_dsp = 1;
                    centroid->phases = {
                        {{ 0,  0}, ndiv_dsp},
                        {{ 0,  0}, ndiv_dsp},
                        {{ 0,  0}, ndiv_dsp},
                        {{-1, -1}, 1},
                        {{ 1,  1}, ndiv_dsp},
                        {{ 1,  1}, ndiv_dsp},
                        {{ 1,  1}, ndiv_dsp},
                        {{ 1,  1}, 1}
                    };
                }
            }
            
            goalTime = 1.8;
            vec3_t jumpPos(0.5, 0.0, comHeight + 0.2);
            vec3_t jumpOri(0.0, deg_to_rad(180.0), 0.0);
            int N = centroid->phases.size()-1;
            dt = goalTime/((real_t)N);
            centroid->waypoints.resize(N+1);
            {
                Centroid::Waypoint& wp = centroid->waypoints[0];
                wp.value  = Centroid::Waypoint::Value (0.0, dt, startPos, startOri, zero3, zero3);
                wp.weight = Centroid::Waypoint::Weight(1.0, 1, 10.0*one3, 10.0*one3, 10.0*one3, 10.0*one3);
                wp.ends.resize(nend);
                for(int iend = 0; iend < nend; iend++){
		            wp.ends[iend].value  = Centroid::Waypoint::End::Value (startPos + endConf[iend].basePos + endConf[iend].posOrigin, zero3, zero3, zero3, infi);
                    wp.ends[iend].weight = Centroid::Waypoint::End::Weight(10.0*one3, 10.0*one3, 10.0*one3, 10.0*one3, 1.0, 1.0*one2, 1.0*one3);
                }
            }
            {
                Centroid::Waypoint& wp = centroid->waypoints[1];
                wp.value  = Centroid::Waypoint::Value (inf, inf, startPos, startOri, inf3, zero3);
                wp.weight = Centroid::Waypoint::Weight(1, 1, 1.0*one3, one3, 1.0*one3, one3);
                wp.ends.resize(nend);
                for(int iend = 0; iend < nend; iend++){
		            wp.ends[iend].value  = Centroid::Waypoint::End::Value (startPos + endConf[iend].basePos + endConf[iend].posOrigin, zero3, zero3, zero3, infi);
                    wp.ends[iend].weight = Centroid::Waypoint::End::Weight(1*one3,1*one3,1*one3,1*one3, 1.0, 1*one2, 1*one3);
                }
            }
            {
                Centroid::Waypoint& wp = centroid->waypoints[2];
                wp.value  = Centroid::Waypoint::Value (inf, inf, startPos, startOri, inf3, zero3);
                wp.weight = Centroid::Waypoint::Weight(inf, inf, 1.0*one3, one3, 1.0*one3, one3);
                wp.ends.resize(nend);
                for(int iend = 0; iend < nend; iend++){
		            wp.ends[iend].value  = Centroid::Waypoint::End::Value (startPos + endConf[iend].basePos + endConf[iend].posOrigin, zero3, zero3, zero3, infi);
                    wp.ends[iend].weight = Centroid::Waypoint::End::Weight(1*one3,1*one3,1*one3,1*one3, 1.0, 1*one2, 1*one3);
                }
            }
            {
                Centroid::Waypoint& wp = centroid->waypoints[3];
                wp.value  = Centroid::Waypoint::Value (inf, inf, startPos, startOri - vec3_t(0.0, deg_to_rad(45.0), 0.0), zero3, vec3_t(0.0, -10.0, 0.0));
                wp.weight = Centroid::Waypoint::Weight(inf, inf, one3, one3, one3, one3);
                wp.ends.resize(nend);
                for(int iend = 0; iend < nend; iend++){
		            wp.ends[iend].value  = Centroid::Waypoint::End::Value (startPos + endConf[iend].basePos + endConf[iend].posOrigin, zero3, zero3, zero3, infi);
                    wp.ends[iend].weight = Centroid::Waypoint::End::Weight(one3, one3, one3, one3, 1.0, 0.1*one2, one3);
                }
            }
            {
                Centroid::Waypoint& wp = centroid->waypoints[4];
                wp.value  = Centroid::Waypoint::Value (inf, inf, goalPos, goalOri + vec3_t(0.0, deg_to_rad(45.0), 0.0), zero3, vec3_t(0.0, -10.0, 0.0));
                wp.weight = Centroid::Waypoint::Weight(inf, inf, one3, one3, one3, one3);
                wp.ends.resize(nend);
                for(int iend = 0; iend < nend; iend++){
		            wp.ends[iend].value  = Centroid::Waypoint::End::Value (goalPos + endConf[iend].basePos + endConf[iend].posOrigin, zero3, zero3, zero3, infi);
                    wp.ends[iend].weight = Centroid::Waypoint::End::Weight(one3, one3, one3, one3, 1.0, one2, one3);
                }
            }
            {
                Centroid::Waypoint& wp = centroid->waypoints[5];
                wp.ends.resize(nend);
                wp.value  = Centroid::Waypoint::Value (inf, inf, goalPos, goalOri, inf3, zero3);
                wp.weight = Centroid::Waypoint::Weight(inf, inf, one3, one3, one3, one3);
                for(int iend = 0; iend < nend; iend++){
		            wp.ends[iend].value  = Centroid::Waypoint::End::Value (goalPos + endConf[iend].basePos + endConf[iend].posOrigin, zero3, zero3, zero3, infi);
                    wp.ends[iend].weight = Centroid::Waypoint::End::Weight(1*one3,1*one3,1*one3,1*one3, 1.0, 1*one2, 1*one3);
                }
            }
            {
                Centroid::Waypoint& wp = centroid->waypoints[6];
                wp.value  = Centroid::Waypoint::Value (dt*N, dt, goalPos, goalOri, zero3, zero3);
                wp.weight = Centroid::Waypoint::Weight(1, 1, 10.0*one3, 10.0*one3, 10.0*one3, 10.0*one3);
                wp.ends.resize(nend);
                for(int iend = 0; iend < nend; iend++){
		            wp.ends[iend].value  = Centroid::Waypoint::End::Value (goalPos + endConf[iend].basePos + endConf[iend].posOrigin, zero3, zero3, zero3, infi);
                    wp.ends[iend].weight = Centroid::Waypoint::End::Weight(10.0*one3, 10.0*one3, 10.0*one3, 10.0*one3, 1.0, 1.0*one2, 1.0*one3);
                }
            }
        }
        centroid->ends.resize(nend);
	    for(int iend = 0; iend < nend; iend++){
            centroid->ends[iend].basePos = endConf[iend].basePos;
            centroid->ends[iend].posMin  = endConf[iend].posMin;
            centroid->ends[iend].posMax  = endConf[iend].posMax; 
		    centroid->ends[iend].copMin  = vec2_t(-0.05, -0.01);
		    centroid->ends[iend].copMax  = vec2_t( 0.05,  0.01);

            centroid->ends[iend].stiffnessMax = endConf[iend].stiffMax;
            centroid->ends[iend].cmpOffset  = endConf[iend].cmpOffset;
            centroid->ends[iend].lockOri    = false;
            centroid->ends[iend].lockCmp    = false;
            centroid->ends[iend].lockMoment = false;
	    }
                
        int N = centroid->NumSteps();
	    for(int k = 0; k <= N; k++)
		    new Tick(world.get(), k, k*dt);
		
        centroid->SetScaling();
	    world->Init();

        centroid->Setup();
        centroid->Reset(true, true, true);
        centroid->Prepare();
        
	    //graph->solver->Enable(ID(DiMP::ConTag::CentroidPosT      ), false);
	    //graph->solver->Enable(ID(DiMP::ConTag::CentroidPosR      ), false);
	    //graph->solver->Enable(ID(DiMP::ConTag::CentroidVelT      ), false);
	    //graph->solver->Enable(ID(DiMP::ConTag::CentroidVelR      ), false);
	    //graph->solver->Enable(ID(DiMP::ConTag::CentroidTime      ), false);
        //graph->solver->Enable(ID(DiMP::ConTag::CentroidEndPos    ), false);
	    //graph->solver->Enable(ID(DiMP::ConTag::CentroidEndVel    ), false);
	    //graph->solver->Enable(ID(DiMP::ConTag::CentroidEndStiff  ), false);
	    world->solver->Enable(ID(ConTag::CentroidEndPosRange  ), false);
	    world->solver->Enable(ID(ConTag::CentroidEndContact), false);
        world->solver->Enable(ID(ConTag::CentroidEndFriction), false);
	    world->solver->Enable(ID(ConTag::CentroidEndMomentRange), false);
        
        world->solver->SetCorrection(ID(), 0.5);
        world->solver->param.method = Solver::Method::DDP;
	    world->solver->param.regularization = 10;
	    world->solver->param.stateRegularization = 10;
        world->solver->param.hastyStepSize  = false;
	    world->solver->param.cutoffStepSize = 0.1;
	    world->solver->param.minStepSize    = 1.0;
	    world->solver->param.maxStepSize    = 1.0;
        world->solver->param.verbose        = true;
        world->solver->param.parallelize    = false;
        world->solver->param.fixInitialState = true;
    
        fileDuration = fopen("duration.csv", "w");
        fileCost     = fopen("cost.csv", "w");

        fprintf(fileCost,
            "iter, step, obj, Tpre, Tdir, Tstep, Tmod, Trans, Tcost, Tcostgrad, Tback, Tpre, Tprestep, Tstep, Tfinish, wcompl\n"
        );
    }

    virtual void OnStep(){
        fprintf(fileCost,
            "%d, %f, %f, %d, %d, %d, %d, "
            "%d, %d, %d, %d, "
            "%d, %d, %d, %d, "
            "\n",
            world->solver->status.iterCount,
		    world->solver->status.stepSize,
		    world->solver->status.obj,
		    world->solver->status.timePre,
		    world->solver->status.timeDir,
		    world->solver->status.timeStep,
		    world->solver->status.timeMod,
            world->solver->status.timeTrans,
            world->solver->status.timeCost,
            world->solver->status.timeCostGrad,
            world->solver->status.timeBack,
            world->TPrepare,
		    world->TPrepareStep,
		    world->TStep,
		    world->TFinish);
       fflush(fileCost);
    }

    void SavePlan(){
        static int idx = 0;
        char filename[256];
        sprintf(filename, "plan_centroid.csv");
        FILE* file = fopen(filename, "w");
        idx++;

        fprintf(file, 
            "k, "
            "time, duration, "
            "cen_pos_t_x, cen_pos_t_y, cen_pos_t_z, "
            "cen_vel_t_x, cen_vel_t_y, cen_vel_t_z, "
            "cen_pos_r_x, cen_pos_r_y, cen_pos_r_z, "
            "cen_vel_r_x, cen_vel_r_y, cen_vel_r_z, "
        );
        for(int i = 0; i < centroid->ends.size(); i++){
            fprintf(file,
                "end%d_pos_t_x, end%d_pos_t_y, end%d_pos_t_z, "
                "end%d_vel_t_x, end%d_vel_t_y, end%d_vel_t_z, "
                "end%d_pos_r_x, end%d_pos_r_y, end%d_pos_r_z, "
                "end%d_vel_r_x, end%d_vel_r_y, end%d_vel_r_z, "
                "end%d_force_t_x, end%d_force_t_y, end%d_force_t_z, "
                "end%d_force_r_x, end%d_force_r_y, end%d_force_r_z, ",
                i, i, i,
                i, i, i,
                i, i, i,
                i, i, i,
                i, i, i,
                i, i, i
                );
        }
        for(int i = 0; i < centroid->ends.size(); i++){
            fprintf(file,
                "end%d_iface, ",
                i
                );
        }
        fprintf(file, "\n");

        for(int k = 0; k < world->ticks.size(); k++){
            auto  key = (CentroidKey*)centroid->traj.GetKeypoint(k);
            auto& d   = key->data;
        
            fprintf(file,
                "%d, "
                "%f, %f, "
                "%f, %f, %f, "
                "%f, %f, %f, ",
                k, 
                d.time, d.duration,
                d.pos_t.x(), d.pos_t.y(), d.pos_t.z(), 
                d.vel_t.x(), d.vel_t.y(), d.vel_t.z()
            );
            fprintf(file,
                "%f, %f, %f, "
                "%f, %f, %f, ",
                d.pos_r.x(), d.pos_r.y(), d.pos_r.z(), 
                d.L.x(), d.L.y(), d.L.z()
            );
            for(int i = 0; i < key->ends.size(); i++){
                auto& dend = d.ends[i];

                fprintf(file,
                    "%f, %f, %f, "
                    "%f, %f, %f, ",
                    dend.pos_t.x(), dend.pos_t.y(), dend.pos_t.z(), 
                    dend.vel_t.x(), dend.vel_t.y(), dend.vel_t.z()
                );
                fprintf(file,
                    "%f, %f, %f, "
                    "%f, %f, %f, ",
                    dend.pos_r.x(), dend.pos_r.y(), dend.pos_r.z(), 
                    dend.vel_r.x(), dend.vel_r.y(), dend.vel_r.z()
                );
                fprintf(file,
                    "%f, %f, %f, "
                    "%f, %f, %f, ",
                    dend.force_t.x(), dend.force_t.y(), dend.force_t.z(), 
                    dend.force_r.x(), dend.force_r.y(), dend.force_r.z()
                );
            }
            for(int i = 0; i < key->ends.size(); i++){
                fprintf(file,
                    "%d, ",
                    key->data_des.ends[i].iface
                );    
            }
            fprintf(file, "\n");
        }

        fclose(file);
    }

    void SaveTraj(){
        char filename[256];
        sprintf(filename, "traj_centroid.csv");
        FILE* file = fopen(filename, "w");
    
        fprintf(file, 
            "k, "
            "cen_pos_t_x, cen_pos_t_y, cen_pos_t_z, "
            "cen_pos_r_w, cen_pos_r_x, cen_pos_r_y, cen_pos_r_z, "
            "cen_vel_t_x, cen_vel_t_y, cen_vel_t_z, "
            //"cen_vel_r_x, cen_vel_r_y, cen_vel_r_z, "
            "cen_L_x, cen_L_y, cen_L_z, "
        );
        for(int i = 0; i < centroid->ends.size(); i++){
            fprintf(file,
                "end%d_pos_t_x, end%d_pos_t_y, end%d_pos_t_z, "
                "end%d_pos_r_w, end%d_pos_r_x, end%d_pos_r_y, end%d_pos_r_z, "
                "end%d_vel_t_x, end%d_vel_t_y, end%d_vel_t_z, "
                "end%d_vel_r_x, end%d_vel_r_y, end%d_vel_r_z, "
                "end%d_force_t_x, end%d_force_t_y, end%d_force_t_z, "
                "end%d_force_r_x, end%d_force_r_y, end%d_force_r_z, ",
                i, i, i,
                i, i, i, i,
                i, i, i,
                i, i, i,
                i, i, i,
                i, i, i
                );
        }
        fprintf(file, "\n");


        auto key = (CentroidKey*)centroid->traj.GetKeypoint(world->ticks.back()->idx);
        real_t tf = key->var_time->val;
        const real_t dt = 0.01;
        CentroidData d;

        for(real_t t = 0.0; t <= tf; t += dt){
            centroid->CalcState(t, d);

            fprintf(file,
                "%f, "
                "%f, %f, %f, "
                "%f, %f, %f, %f, "
                "%f, %f, %f, "
                //"%f, %f, %f, "
                "%f, %f, %f, ",
                t, 
                d.pos_t.x(), d.pos_t.y(), d.pos_t.z(), 
                d.pos_r.w(), d.pos_r.x(), d.pos_r.y(), d.pos_r.z(),
                d.vel_t.x(), d.vel_t.y(), d.vel_t.z(), 
                d.L.x(), d.L.y(), d.L.z()
            );
            //vec3_t mom;
            for(int i = 0; i < key->ends.size(); i++){
                //mom += (d.ends[i].pos_t - d.pos_t) % d.ends[i].force_t + d.ends[i].force_r;
            
                fprintf(file,
                    "%f, %f, %f, "
                    "%f, %f, %f, %f, "
                    "%f, %f, %f, "
                    "%f, %f, %f, "
                    "%f, %f, %f, "
                    "%f, %f, %f, ",
                    d.ends[i].pos_t.x(), d.ends[i].pos_t.y(), d.ends[i].pos_t.z(), 
                    d.ends[i].pos_r.w(), d.ends[i].pos_r.x(), d.ends[i].pos_r.y(), d.ends[i].pos_r.z(),
                    d.ends[i].vel_t.x(), d.ends[i].vel_t.y(), d.ends[i].vel_t.z(), 
                    d.ends[i].vel_r.x(), d.ends[i].vel_r.y(), d.ends[i].vel_r.z(),
                    d.ends[i].force_t.x(), d.ends[i].force_t.y(), d.ends[i].force_t.z(),
                    d.ends[i].force_r.x(), d.ends[i].force_r.y(), d.ends[i].force_r.z()
                );    
            }
            //fprintf(file, "%f, %f, %f, ", mom.x, mom.y, mom.z);
            fprintf(file, "\n");
        }

        fclose(file);
    }

};

}

// main function
int main(int argc, const char** argv) {
    dymp::MyApp app;
    app.Init();
    app.Loop();
    app.Cleanup();

    return 1;
}
