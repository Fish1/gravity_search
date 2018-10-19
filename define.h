#define DIM 3
// number of objects around the sun
#define POPULATION_SIZE 25
#define DOMAIN 8.0
// mass of the sun, one magnitude higher than the largest particle seems to work
#define SUN_MASS 40000.0
// if a particle is at the current best fitness it will be this size
#define MAX_PARTICLE_MASS 2000.0
// if a particle is at the 
#define MIN_PARTICLE_MASS 200.0
// how far the sun can wiggle
#define WIGGLE_STEP 0.01
// this currently doesn't do what it says...
#define MAX_MOVEMENT 1.0
// kill a portion of the population after this many fitness calls
#define FITNESS_CALLS_TO_KILL 5000
// durring a kill, kill this many particles
#define NUMBER_TO_KILL 5

// how many vertices to put in between a unit line
#define PER_UNIT 3
// render a frame after this many update cylces
#define UPDATES_PER_FRAME 1
