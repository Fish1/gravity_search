#include "ofApp.h"

#include <vector>

#include <iostream>

// Fitness function
double ofApp::function(const double * coords, unsigned int dim)
{
	++m_fitnessCalls;
	++m_eraseCounter;

	double sum1 = 0;

	double sum2 = 0;

	for(unsigned int index = 1; index <= dim; ++index)
	{
		sum1 += pow((coords[index - 1] + pow(-1.0, index) * (index % 4)), 2);
		
		sum2 += pow(coords[index - 1], index);	
	}

	double result = (-1.0) * sqrt(sum1) + sin(sum2);

	return result;
}

/*
double ofApp::function(const double * coords, unsigned int dim)
{
	++m_fitnessCalls;
	++m_eraseCounter;

	double sum1 = 0;

	for(unsigned int index = 1; index <= dim; ++index)
	{
		sum1 += pow(coords[index - 1], 2);
	}

	double result = sin(sqrt(sum1)) + sqrt(abs(coords[0]));

	return result;
}*/

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
			double currentY = function(coord, 2);
			
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

void ofApp::createMesh3D(unsigned int exp)
{
	mesh.clear();

	double bestFitness = -DBL_MAX;
	double worstFitness = DBL_MAX;

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

		for(int y = 0; y < checks; ++y)
		{
			double currentY = ((double)y / (double)perUnit) - ((double)size / 2.0);

			for(int x = 0; x < checks; ++x)
			{
				// the x position of the current vertex
				double currentX = ((double)x / (double)perUnit) - ((double)size / 2.0);

				// pass in these coordinates to the fitness function to get the y position
				double coord [] = {currentX, currentY, currentZ};

				// the y position of the current vertex
				double fitness = function(coord, 3);
			
				if(fitness > bestFitness)
					bestFitness = fitness;

				if(fitness < worstFitness)
					worstFitness = fitness;
	
				ofVec3f point(currentX, currentY, currentZ);
				mesh.addVertex(point);
			}
		}
	}

	bestFitness = bestFitness + abs(worstFitness);

	// Create Colors
	for(int z = 0; z < checks; ++z)
	{
		// the z position of the current vertex
		double currentZ = ((double)z / (double)perUnit) - ((double)size / 2.0);

		for(int y = 0; y < checks; ++y)
		{
			double currentY = ((double)y / (double)perUnit) - ((double)size / 2.0);

			for(int x = 0; x < checks; ++x)
			{
				// the x position of the current vertex
				double currentX = ((double)x / (double)perUnit) - ((double)size / 2.0);

				// pass in these coordinates to the fitness function to get the y position
				double coord [] = {currentX, currentY, currentZ};

				// the y position of the current vertex
				double fitness = function(coord, 3);
	
				fitness += abs(worstFitness);
				
				double closeness = fitness / bestFitness;
	
				closeness = pow(closeness, exp);	
		
				// 255 255 255   0
				// 255 255 255 255
	
				double r =   0 * closeness + 255;
				double g =   0 * closeness + 255;
				double b =   0 * closeness + 255;
				double a = 255 * closeness + 0;	
				
				mesh.addColor(ofColor(r,g,b,a));
			}
		}
	}
}

void ofApp::createMesh3D_test(unsigned int exp)
{
	mesh.clear();

	double bestFitness = -DBL_MAX;
	double worstFitness = DBL_MAX;

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
		double coord [] = {x, y, z};

		// the y position of the current vertex
		double fitness = function(coord, 3);
	
		if(fitness > bestFitness)
			bestFitness = fitness;

		if(fitness < worstFitness)
			worstFitness = fitness;

		ofVec3f point(x, y, z);
		mesh.addVertex(point);
	}

	bestFitness += worstFitness;

	auto vertices = mesh.getVertices();

	for(auto & vert : vertices)
	{
		double coord [] = {vert.x, vert.y, vert.z};

		double fitness = function(coord, 3);

		fitness += worstFitness;

		double closeness = fitness / bestFitness;
	
		if(closeness >= 0.0)	
			closeness = pow(closeness, exp);	
		else if(closeness < 0.0)
			closeness = 0.0;	

		// 255   0   0 100
		// 255 255 255 255

		double r =   0.0 * closeness + 255.0;
		double g = 255.0 * closeness + 0.0;
		double b = 255.0 * closeness + 0.0;
		double a = 255.0 * closeness + 0.0;	
		
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
		createMesh3D_test(m_exp);
	}


	// Initialize the camera closer to our graph
	cam.setTarget(glm::vec3(0.0f,-5.0f,0.0f));
	cam.setDistance(20.0f);

	// Initialize the population

	for(unsigned int i = 0; i < POPULATION_SIZE; ++i)
	{
		double startPosition[DIM];
		for(unsigned int d = 0; d < DIM; ++d)
			startPosition[d] = domain(gen);

		double startVelocity[DIM];
		for(unsigned int d = 0; d < DIM; ++d)
			startVelocity[d] = 0;

		Particle p(DIM, startPosition, startVelocity, 1000.0);

		m_particles.push_back(p);
	}
	
	double startPosition[DIM];
	for(unsigned int d = 0; d < DIM; ++d)
		startPosition[d] = domain(gen);
	
	double startVelocity[DIM];
	for(unsigned int d = 0; d < DIM; ++d)
		startVelocity[d] = 0;

	m_sun = Particle(DIM, startPosition, startVelocity, SUN_MASS);
	
	for(unsigned int d = 0; d < DIM; ++d)
		std::cout << m_sun.getPosition(d) << " ";

	std::cout << std::endl;

	// Initialize the clock
	m_start = std::clock();

	m_bestFitness = -DBL_MAX;
}

//--------------------------------------------------------------
void ofApp::update()
{
	//if(((std::clock() - m_start) / (double) CLOCKS_PER_SEC) < 0.0000001)
	{
	//	return;
	}

	++m_updateCount;
/*
	for(unsigned int d = 0; d < DIM; ++d)
		std::cout << m_sun.getPosition(d) << " ";

	std::cout << std::endl;*/

	if(m_eraseCounter > 5000)
	{
		m_particles.erase(m_particles.begin(), m_particles.begin() + 5);
		
		for(unsigned int i = 0; i < 5; ++i)
		{
			double startPosition[DIM];
			for(unsigned int d = 0; d < DIM; ++d)
				startPosition[d] = domain(gen);

			double startVelocity[DIM];
			for(unsigned int d = 0; d < DIM; ++d)
				startVelocity[d] = 0;

			Particle p(DIM, startPosition, startVelocity, 1000.0);

			m_particles.push_back(p);
		}

		m_eraseCounter = 0;
	}

	m_sun.wiggle();

	if(m_sun.isOutOfBounds())
		m_sun.revertPosition();

	double fitness = function(m_sun.getPosition(), DIM);

	m_sun.setCurrentFitness(fitness);

	for(Particle & p : m_particles)
	{
		p.updateVelocity(m_sun);

		for(Particle & pa : m_particles)
		{
			if(&p == &pa)
				continue;
				
			p.updateVelocity(pa);
		}

		p.updatePosition();

		if(p.isOutOfBounds())
			continue;
		
		fitness = function(p.getPosition(), DIM);

		p.setCurrentFitness(fitness);

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
		}
		else if(fitness > m_sun.getCurrentFitness())
		{
			m_sun.setPosition(p.getPosition());
			
			m_sun.setCurrentFitness(p.getCurrentFitness());
		}
		
		double fitnessDiff = m_bestFitness - fitness;
		fitnessDiff = abs(fitnessDiff);
		p.setMass(10000.0 - 
			(9000.0 * (fitnessDiff / (fitnessDiff + 1))));
	}

	m_start = std::clock();
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
void ofApp::keyPressed(int key){

	if(DIM == 3)
	{
		if(key == OF_KEY_RIGHT)
			m_exp += 1;
		else if(key == OF_KEY_LEFT)
		{
			m_exp -= 1;
			if(m_exp < 1)
				m_exp = 1;	
		}
			
		createMesh3D_test(m_exp);
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