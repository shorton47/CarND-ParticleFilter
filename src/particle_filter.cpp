/*
 * particle_filter.cpp
 *
 *  Originally Created on: Dec 12, 2016, Author: Tiffany Huang
 *  Filled in: v1   6/24/17     SWH
 *
 */
#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

// Class global for debug
bool PDEBUG = true;




//---------------
// For each iteration of the master Particle Filter engine (from main.cpp), there are up to 4 main methods called:
//
//   1. init
//   2. predicition
//   3. updateWeights
//   4. resample
//
//   There is also a Helper Methods section below that is divided into the following parts:
//
//---------------



// Section #1
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    /**
     * init Initializes particle filter by initializing particles to Gaussian
     *   distribution around first position and all the weights to 1.
     * @param x Initial x position [m] (simulated estimate from GPS)
     * @param y Initial y position [m] (simulated estimate from GPS)
     * @param theta Initial orientation [rad]
     * @param std[] Array of dimension 3 [standard deviation of x [m], std of y [m], std of yaw [rad]]
     */

    Particle one_particle;      // single particle object
    default_random_engine gen;  // rv sequence object
    
    
    // Set # of particles
    num_particles = 5;
    if (PDEBUG) cout << "Number of particle to be used=" << num_particles << endl;
    
    // Normal (Gaussian) distributions with means position and yaw w/ st deviations's of std[].
    normal_distribution <double> dist_x    (x,     std[0]);
    normal_distribution <double> dist_y    (y,     std[1]);
    normal_distribution <double> dist_theta(theta, std[2]);
    
    // Generate particles w/ position with additive gaussian noise. "gen" is the random engine initialized earlier
    for (int i = 0; i < num_particles; ++i) {
        
        one_particle.id     = i+1;
        one_particle.x      = dist_x(gen);
        one_particle.y      = dist_y(gen);
        one_particle.theta  = dist_theta(gen);
        one_particle.weight = 1.0F;
        
        // Add this particle to array of particles.
        particles.push_back(one_particle);
        weights.push_back(1.0F);  // <TODO> What is diff between one_particle.weight and weights?????
        
        // Print particle
        if (PDEBUG) cout << "Particle: " << one_particle.id << " " << one_particle.x     << " "
                                         << one_particle.y  << " " << one_particle.theta << " " <<  endl;
        }
    
    is_initialized = true;
    
} // init


// Section #2
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


// 3 parts to predicition
    
    //P-1 Sample
    //P-2 Calculate weight
    //P-3 Update state vector



}


// Section #3
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html
}


// Section #4
void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
}



//---------------
// Helper Methods
//
//
//---------------

// Called by updateWeights
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}




Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}


string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
