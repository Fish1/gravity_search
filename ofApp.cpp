#include "ofApp.h"

#include <vector>

#include <iostream>

// Fitness function
/*
double ofApp::function(const double * coords, unsigned int dim, bool count = true)
{
	if(count)
	{
		++m_fitnessCalls;
		++m_eraseCounter;
	}

	double sum1 = 0;

	double sum2 = 0;

	for(unsigned int index = 1; index <= dim; ++index)
	{
		sum1 += pow((coords[index - 1] + pow(-1.0, index) * (index % 4)), 2);
		
		sum2 += pow(coords[index - 1], index);	
	}

	double result = (-1.0) * sqrt(sum1) + sin(sum2);

	return result;
}*/

double ofApp::function(const double * coords, unsigned int dim, bool count = true)
{
	if(count)
	{
		++m_fitnessCalls;
		++m_eraseCounter;
	}

	double distance = 0;
	double tmpCoord[dim];
	for(unsigned int d = 0; d < dim; ++d)
	{
		tmpCoord[d] = coords[d] - 1.0;
	}

	for(unsigned int d = 0; d < dim; ++d)
	{
		distance += pow(tmpCoord[d], 2.0);
	}

	distance = sqrt(distance);
	
	double underCos = 8 * distance;

	double numerator = cos(underCos);
	double denominator = distance + 0.1;

	return numerator / denominator;
}

void ofApp::createMesh2D()
{
	// size is from -8 to 8
	const int size = DOMAIN * 2;
	// how many vertices per 1 unit
	const int perUnit = PER_UNIT;
	// square root of the number of vertices
	const int checks = perUnit * size;

	// Create Verticies
	for(int z = 0; z < checks; ++z)
	{
		// the z position of the current vertex
		double currentZ = ((double)z / (double)perUnit) - ((double)size / 2.0);

		for(int x = 0; x < checks; ++x)
		{
			// the x position of the current vertex
			double currentX = ((double)x / (double)perUnit) - ((double)size / 2.0);

			// pass in these coordinates to the fitness function to get the y position
			double coord [] = {currentX, currentZ};

			// the y position of the current vertex
			double currentY = function(coord, 2, false);
			
			ofVec3f point(currentX, currentY, currentZ);
			mesh.addVertex(point);
		}
	}

	// Create indices

	for(unsigned int y = 0; y < checks - 1; ++y)
	{
		for(unsigned int x = 0; x < checks; ++x)
		{
			unsigned int current = x + checks * y;
			unsigned int below = x + checks * (y + 1);
			unsigned int left = (x - 1) + checks * y;
			unsigned int belowRight = (x + 1) + checks * (y + 1);

			if(x == 0)
			{
				mesh.addIndex(current);
				mesh.addIndex(below);
				mesh.addIndex(belowRight);	
			}
			else if(x == checks - 1)
			{
				mesh.addIndex(current);
				mesh.addIndex(left);
				mesh.addIndex(below);
			}
			else
			{
				mesh.addIndex(current);
				mesh.addIndex(below);
				mesh.addIndex(belowRight);
				
				mesh.addIndex(current);
				mesh.addIndex(left);
				mesh.addIndex(below);
			}
		}
	}
}

void ofApp::createMesh3D()
{
	mesh.clear();

	// size is from -8 to 8
	const int size = DOMAIN * 2;
	// how many vertices per 1 unit
	const int perUnit = PER_UNIT;
	// square root of the number of vertices
	unsigned int checks = perUnit * size;

	checks = pow(checks, 3);
	

	for(unsigned int index = 0; index < checks; ++index)
	{
		double x = domain(gen);
		double y = domain(gen);
		double z = domain(gen);
		
		// pass in these coordinates to the fitness function to get the y position
		// double coord [] = {x, y, z};

		// the y position of the current vertex
		// double fitness = function(coord, 3);
	
		ofVec3f point(x, y, z);
		mesh.addVertex(point);
	}
}

void ofApp::colorMesh3D(unsigned int exp)
{
	mesh.clearColors();

	auto vertices = mesh.getVertices();

	for(auto & vert : vertices)
	{
		double coord [] = {vert.x, vert.y, vert.z};

		double fitness = function(coord, 3, false);

		// makes it 0 or above
		if(fitness < m_worstFitness)
			m_worstFitness = fitness;

		fitness += abs(m_worstFitness);

		double tmpBestFitness = m_bestFitness;
		// because we shifted the fitness up
		// we need to shift the best fitness up
		tmpBestFitness += abs(m_worstFitness);

		double closeness = fitness / tmpBestFitness;

		closeness = pow(closeness, (double)(exp));

		if(exp == 1)
		{
			closeness = pow(closeness, 2.0);
		}

		// 255   0   0   0
		// 255 255 255 255

		double r =   0.0 * closeness + 255.0;
		double g = 255.0 * closeness + 0.0;
		double b = 255.0 * closeness + 0.0;
		double a = 255.0 * closeness + 0.0;	
		if(exp == 1)
			a = 255.0;
		
		mesh.addColor(ofColor(r,g,b,a));
	}

}

//--------------------------------------------------------------
void ofApp::setup()
{
	ofEnableDepthTest();
	ofEnableAlphaBlending();
	ofEnableBlendMode(OF_BLENDMODE_ALPHA);
	// only reset the background when i tell it too
	ofSetBackgroundAuto(false);

	domain = std::uniform_real_distribution<double>(-DOMAIN, DOMAIN);

	if(DIM == 2)
	{
		createMesh2D();
	}
	else if(DIM == 3)
	{
		m_exp = 2;
		createMesh3D();
	}


	// Initialize the camera closer to our graph
	cam.setTarget(glm::vec3(0.0f,-5.0f,0.0f));
	cam.setDistance(20.0f);

	// Initialize the population
	double startPosition[DIM];
	double startVelocity[DIM];

	for(unsigned int i = 0; i < POPULATION_SIZE; ++i)
	{
		for(unsigned int d = 0; d < DIM; ++d)
			startPosition[d] = domain(gen);

		for(unsigned int d = 0; d < DIM; ++d)
			startVelocity[d] = 0;

		Particle p(DIM, startPosition, startVelocity, 1000.0);

		m_particles.push_back(p);
	}
	
	for(unsigned int d = 0; d < DIM; ++d)
		startPosition[d] = domain(gen);
	
	for(unsigned int d = 0; d < DIM; ++d)
		startVelocity[d] = 0;

	m_sun = Particle(DIM, startPosition, startVelocity, SUN_MASS);
	
	std::cout << std::endl;


	// Initialize worst and best fitness
	m_worstFitness = DBL_MAX;	
	m_bestFitness = -DBL_MAX;
}

//--------------------------------------------------------------
void ofApp::update()
{
	++m_updateCount;

	// Reset some of the particles after a while
	if(m_eraseCounter > FITNESS_CALLS_TO_KILL)
	{
		m_particles.erase(m_particles.begin(), m_particles.begin() + NUMBER_TO_KILL);
		
		double startPosition[DIM];
		double startVelocity[DIM];

		for(unsigned int i = 0; i < NUMBER_TO_KILL; ++i)
		{
			for(unsigned int d = 0; d < DIM; ++d)
				startPosition[d] = domain(gen);

			for(unsigned int d = 0; d < DIM; ++d)
				startVelocity[d] = 0;

			Particle p(DIM, startPosition, startVelocity, 1000.0);

			m_particles.push_back(p);
		}

		m_eraseCounter = 0;
	}

	// wiggle the sun a bit
	m_sun.wiggle();
	if(m_sun.isOutOfBounds())
		m_sun.revertPosition();
	double fitness = function(m_sun.getPosition(), DIM);
	m_sun.setCurrentFitness(fitness);

	// update each particle
	for(Particle & p : m_particles)
	{
		// the particle reacts to the suns gravity
		p.updateVelocity(m_sun);

		for(Particle & pa : m_particles)
		{
			if(&p == &pa)
				continue;
			
			// the particle reacts to other particles gravity	
			p.updateVelocity(pa);
		}

		const double * prePos = p.getPosition();
		p.updatePosition();
		// if the particle has moved out of bounds, forget about it
		if(p.isOutOfBounds())
		{
			continue;
		}
		const double * curPos = p.getPosition();

		// calculate the distance the particle moved
		double distanceMoved = 0;	
		for(unsigned int d = 0; d < DIM; ++d)
		{
			distanceMoved += pow(curPos[d] - prePos[d], 2);
		}
		distanceMoved = sqrt(distanceMoved);

		// get the particles fitness 	
		fitness = function(p.getPosition(), DIM);
		p.setCurrentFitness(fitness);

		// if the fitness is greater than the best fitness
		// update the best fitness and move the sun there
		if(fitness > m_bestFitness)
		{
			m_bestFitness = fitness;

			for(unsigned int d = 0; d < DIM; ++d)
				m_bestPosition[d] = p.getPosition(d);

			m_sun.setPosition(m_bestPosition);
			
			m_sun.setCurrentFitness(m_bestFitness);

			std::cout << "Fitness: " << m_bestFitness << std::endl;
			for(unsigned int d = 0; d < DIM; ++d)
			{
				std::cout << m_bestPosition[d] << " ";
			}
			std::cout << std::endl;
			std::cout << "Calls: " << m_fitnessCalls << std::endl;
			std::cout << std::endl;
	
			// if we are in dimension 3, recolor the mesh	
			if(DIM == 3)
			colorMesh3D(m_exp);
		}
		// if the particles fitness is greater than the suns current fitness
		// move the sun to the particle
		else if(fitness > m_sun.getCurrentFitness())
		{
			m_sun.setPosition(p.getPosition());
			
			m_sun.setCurrentFitness(p.getCurrentFitness());
		}

		// find the worst fitness (for coloring the 3D mesh)
		if(fitness < m_worstFitness)
		{
			m_worstFitness = fitness;
			if(DIM == 3)
			colorMesh3D(m_exp);
		}
	
		// if the particle is at a good fitness increase its mass
		// if the particle is at a bad fitness decrease its mass	
		double fitnessDiff = abs(m_bestFitness - fitness);
		fitnessDiff = abs(fitnessDiff);
		p.setMass(MAX_PARTICLE_MASS - 
			((MAX_PARTICLE_MASS - MIN_PARTICLE_MASS) * (fitnessDiff / (fitnessDiff + 1))));
	}
}

//--------------------------------------------------------------
void ofApp::draw(){

	if(m_updateCount % UPDATES_PER_FRAME != 0)
		return;
	
	cam.begin();

	if(DIM == 2)
	{
		ofBackground(0);
		mesh.enableColors();
		ofSetColor(255,255,255);
		mesh.drawWireframe();
		mesh.disableColors();
	}
	if(DIM == 3)
	{

		ofBackground(0);
		mesh.enableColors();
		ofSetColor(255,255,255);
		mesh.drawVertices();
		mesh.disableColors();
	}

	if(DIM == 2)
	{
		ofSetColor(255,255,0);
		ofDrawSphere(m_sun.getPosition(0), m_sun.getCurrentFitness(), m_sun.getPosition(1), 0.3);

		ofSetColor(255,0,0);
		for(Particle & p : m_particles)
		{
			ofDrawSphere(p.getPosition(0), p.getCurrentFitness(), p.getPosition(1), 0.2);
		}
	}
	else if(DIM == 3)
	{
		ofSetColor(255,255,0);
		ofDrawSphere(m_sun.getPosition(0), m_sun.getPosition(1), m_sun.getPosition(2), 0.3);
		
		ofSetColor(255,0,0);
		for(Particle & p : m_particles)
		{
			ofDrawSphere(p.getPosition(0), p.getPosition(1), p.getPosition(2), 0.2);
		}
	}

	cam.end();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key)
{
	if(DIM == 3)
	{
		if(key == OF_KEY_RIGHT)
		{
			m_exp += 1;
			
			colorMesh3D(m_exp);
		}
		else if(key == OF_KEY_LEFT)
		{
			m_exp -= 1;
			if(m_exp < 1)
				m_exp = 1;	
		
			colorMesh3D(m_exp);
		}
			
	}
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
