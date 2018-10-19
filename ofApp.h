#pragma once

#include "ofMain.h"

#include "particle.h"

#include "define.h"

#include <ctime>
#include <random>

class ofApp : public ofBaseApp
{

	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);

	private:
		ofEasyCam cam;
		ofVboMesh mesh;

		Particle m_sun;
		std::vector<Particle> m_particles;

		double m_bestFitness;

		double m_bestPosition[DIM];

		std::clock_t m_start;

		unsigned int m_fitnessCalls = 0;
		// eraseCounter iterates every time the fitness function is called
		// if the counter is above a threshhold it will erase 5 random
		// particles
		unsigned int m_eraseCounter = 0;
		unsigned int m_updateCount = 0;
	
		std::default_random_engine gen;
		std::uniform_real_distribution<double> domain;	

	private:
		double function(const double * coords, unsigned int dim);

		void createMesh2D();

		unsigned int m_exp = 2;
		void createMesh3D(unsigned int exp);
		void createMesh3D_test(unsigned int exp);
};
