#include "particle.h"

#include <cfloat>

#include <cmath>

#include <random>

#include <iostream>

Particle::Particle() :
	m_bestFitness(-DBL_MAX), dist(-WIGGLE_STEP, WIGGLE_STEP)
{
}

Particle::Particle(unsigned int dim, double * position, double * velocity, double mass) :
	m_mass(mass), m_bestFitness(-DBL_MAX),
	m_currentFitness(-DBL_MAX), dist(-WIGGLE_STEP, WIGGLE_STEP)
{
	for(unsigned int d = 0; d < m_dim; ++d)
	{
		m_position[d] = position[d];
		m_previousPosition[d] = m_position[d];
	}

	for(unsigned int d = 0; d < m_dim; ++d)
		m_velocity[d] = velocity[d];
}

void Particle::setCurrentFitness(double fitness)
{
	m_currentFitness = fitness;

	if(m_currentFitness > m_bestFitness)
	{
		m_bestFitness = m_currentFitness;
	}
}

double Particle::getCurrentFitness()
{
	return m_currentFitness;
}

void Particle::updateVelocity(Particle & sun)
{
	for(unsigned int d = 0; d < m_dim; ++d)
	{
		m_previousPosition[d] = m_position[d];
	}

	// Calculate the force of gravity
	double force = (6.67408 * pow(10, -11)) * sun.getMass() * getMass();
	// should be d^2 not d (get rid of square root)
//	double distance = sqrt(pow(getX() - sun.getX(), 2) + pow(getY() - sun.getY(), 2));
//	if(distance == 0)
		//distance = DBL_EPSILON;

//	std::cout << force << " / " << distance << std::endl;
//	force *= distance;

// 	use this for a constance gravity force between objects
//	force = 0.01;

	// Get the  accleration vector
	double delta[m_dim];
		
	for(unsigned int d = 0; d < m_dim; ++d)
	{
		delta[d] = sun.getPosition(d) - getPosition(d);
	}

	// Normalize the acceleration vector
	double mag = 0.0;

	for(unsigned int d = 0; d < m_dim; ++d)
	{
		mag += (delta[d] * delta[d]);
	}

	mag = sqrt(mag);

	if(mag == 0.0)
		mag = DBL_EPSILON;

	for(unsigned int d = 0; d < m_dim; ++d)
	{
		delta[d] /= mag;
	}

	// Weaken the acceleration vector by a magnitude of 100	
	for(unsigned int d = 0; d < m_dim; ++d)
	{
		delta[d] *= force;
	}

	// Apply the acceleration vector to the velocity vector
	for(unsigned int d = 0; d < m_dim; ++d)
	{
		m_velocity[d] += delta[d];
	}
}

void Particle::updatePosition()
{
	for(unsigned int d = 0; d < m_dim; ++d)
	{
		m_position[d] += m_velocity[d] * 0.1;
	}
}

void Particle::wiggle()
{
	for(unsigned int d = 0; d < m_dim; ++d)
	{
		m_position[d] += dist(gen);
	}
}

void Particle::revertPosition()
{
	for(unsigned int d = 0; d < m_dim; ++d)
	{
		m_position[d] = m_previousPosition[d];
	}
}

const double * const Particle::getPosition()
{
	return m_position;
}

double Particle::getPosition(unsigned int dim)
{
	return m_position[dim];
}

void Particle::setPosition(const double * pos)
{
	for(unsigned int d = 0; d < m_dim; ++d)
	{
		m_position[d] = pos[d];
	}
}

void Particle::setPosition(unsigned int dim, double val)
{
	m_position[dim] = val;
}

double Particle::getMass()
{
	return m_mass;
}

void Particle::setMass(double mass)
{
	m_mass = mass;
}

bool Particle::isOutOfBounds()
{
	for(unsigned int d = 0; d < m_dim; ++d)
	{
		if(m_position[d] > DOMAIN || m_position[d] < -DOMAIN)
		{
			return true;
		}
	}

	return false;
}
