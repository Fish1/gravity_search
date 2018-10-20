#pragma once

#include <random>

#include "define.h"

class Particle
{
private:
	unsigned int m_dim;
	
	double * m_position = nullptr;
	double * m_previousPosition = nullptr;

	double * m_velocity = nullptr;

	double m_mass;

	double m_bestFitness;
	double m_currentFitness;
	
	std::default_random_engine gen;
	std::normal_distribution<double> dist;

private: 
	Particle();

public:
	Particle(unsigned int dim, double * position, double * velocity, double mass, double stdWiggle);
	~Particle();

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
