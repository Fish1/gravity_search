#include "ofMain.h"
#include "ofApp.h"

//========================================================================
int main(int argc, char ** argv)
{
	ofSetupOpenGL(1024,768,OF_WINDOW);			// <-------- setup the GL context

	// this kicks off the running of my app
	// can be OF_WINDOW or OF_FULLSCREEN
	// pass in width and height too:

	std::cout << " -- Function Settings --" << std::endl;

	std::cout << "Dimension: ";

	unsigned int dim;
	
	std::cin >> dim;

	std::cout << " -- Swarm Settings --" << std::endl;

	std::cout << "Swarm Size: ";

	unsigned int swarmSize;

	std::cin >> swarmSize;

	std::cout << "Kill Count: ";
	
	unsigned int killCount;

	std::cin >> killCount;

	std::cout << " -- Sun Settings --" << std::endl;
	
	std::cout << "Standard Deviation Wiggle: ";

	double standardDeviationWiggle;

	std::cin >> standardDeviationWiggle;

	ofRunApp(new ofApp(dim, swarmSize, killCount, standardDeviationWiggle));
}
