#pragma once

#include <random>

#include "define.h"

class Particle
{
private:
	static const unsigned int m_dim = DIM;
	
	double m_position[DIM];
	double m_previousPosition[DIM];

	double m_velocity[DIM];

	double m_mass;

	double m_bestFitness;
	double m_currentFitness;
	
	std::default_random_engine gen;
	std::uniform_real_distribution<double> dist;

public:

	Particle();
	Particle(unsigned int dim, double * position, double * velocity, double mass);

	void setCurrentFitness(double fitness);
	double getCurrentFitness();

	void updateVelocity(Particle & sun);
	void updatePosition();

	void revertPosition();

	void wiggle();

	const double * const getPosition();	
	double getPosition(unsigned int dim);

	void setPosition(const double * pos);
	void setPosition(unsigned int dim, double val);

	const double * const getVelocity();
	
	double getMass();
	
	void setMass(double mass);

	bool isOutOfBounds();
};
