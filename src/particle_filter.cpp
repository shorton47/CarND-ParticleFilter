//---------------
// particle_filter.cpp is the main class that services the particle filter simulation for Udacity project. This is called
// by main.cpp and communicates over Websockets with Udacity simulator to run the Particle Filter on vehicles to
// optimally locate a car given gps, lidar (map observations), and a truth map of landmarks. THe Particle implements
// Bayes conditional probability along with the conditional posterior probability to get an optimal estimate using the
// Bicycle model of vehicle motion
//
// For each iteration of the master engine (from main.cpp), there are up to 4 main methods called in this Particle Filter
// class:
//   1. init              - Initialize particles/weights, etc.. with initial gps coordinates
//   2. predicition       - Predict car location based on position estimate plus heading
//   3. updateWeights     - Update location based on a weighted probability with observed landmarks correlated to a map
//   3A.  dataAssociation - Associates closest landmark w/ each landmark observation from particle/vehicle
//   4. resample          - Update position vectors (particles) of position & heading based on the probability weights
//
//
// Originally Created on: Dec 12,    2016, Author: Tiffany Huang
// Programmed (v1)      : Jun 24-29, 2017, Author: Steve Horton
//---------------

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

//---
// Set class level DEBUG here
bool PDEBUG = false;
//---

//-----
// Step #1 - init
//
// Initialize the particle filter by setting the number of particles, initializing all particles to first position
// (based on estimates of x, y, theta and their uncertainties from GPS) and all weights to 1.
// NOTE: all input is in Map coordinates. initializing particles to Gaussian.
//
// Input:
//  @param x Initial x position [m] (simulated estimate from GPS)
//  @param y Initial y position [m] (simulated estimate from GPS)
//  @param theta Initial orientation [rad]
//  @param std[] Array of dimension 3 [standard deviation of gps x [m], std of y [m], std of yaw [rad]]
//
// Updates:
//  particles - class vector of particles (proposed vehicles)
//  weights - class vector of corresponding particle weights/probabilities
//-----
void ParticleFilter::init(double gps_x, double gps_y, double init_heading, double gps_std[]) {
    
    Particle a_particle;        // temp storage for single particle
    
    
    // Set # of particles (NOTE: CLASS VARIABLE)
    num_particles = 99;
    if (PDEBUG) cout << "PF-Init:Number of particles being generated=" << num_particles << endl;
    
    // Normal (Gaussian) distributions with means position and yaw w/ st deviations's of std[].
    default_random_engine gen;  // random variable sequence object
    normal_distribution <double> dist_x    (gps_x,        gps_std[0]);
    normal_distribution <double> dist_y    (gps_y,        gps_std[1]);
    normal_distribution <double> dist_theta(init_heading, gps_std[2]);
    
    // Generate particles + additive gaussian noise. "gen" is the random # engine initialized earlier
    for (int i = 0; i < num_particles; ++i) {
        
        a_particle.id     = i+1;
        a_particle.x      = dist_x(gen);
        a_particle.y      = dist_y(gen);
        a_particle.theta  = dist_theta(gen);
        a_particle.weight = 1.0F;
        
        particles.push_back(a_particle);  // Add to array of particles
        weights.push_back(1.0F);          // Class array of all particle weights. Needed for particle re-sample step.
        
        //
        if (PDEBUG) cout << "PF-Init:Particle #" << i+1 << " " << a_particle.id << " " << a_particle.x << " "
                                                               << a_particle.y  << " " << a_particle.theta << " " <<  endl;
    }
    is_initialized = true;
    
} // init


//-----
// Step #2 - predicition
//
// Predict vehicle state vector each time step by adding measurements (time, velocity & turning rate) to each particle
// AND some Gaussian noise. STD array are deviations of position and heading angle measurements.
//
// Input:
//  @param delta_t   Time step (sec)
//  @param std_pos[] Array of dimension 3 [standard deviation of pos x [m], std of y [m], std of yaw [rad]]
//  @param velocity  (m/s)
//  @param yaw_rate  (rad/s)
//-----
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	
    double  old_x, new_x, old_y, new_y, old_theta, new_theta;
    double  noise_x, noise_y, noise_theta;
    double  v_div_yawrate;
    bool    yaw_rate_not_zero;
    
    // Create gaussian distributions generators with zero means and measurement stds so can add to prediction
    default_random_engine gen;  // rv sequence object
    normal_distribution <double> dist_x    (0, std_pos[0]);
    normal_distribution <double> dist_y    (0, std_pos[1]);
    normal_distribution <double> dist_theta(0, std_pos[2]);
    
    // Figure out which bicycle model equations will be used for calculation based on yaw rate
    if (abs(yaw_rate) > __DBL_EPSILON__) {
        yaw_rate_not_zero = true;
    } else {
        yaw_rate_not_zero = false;
    }
    
    // Predict all particle's location and heading
    for (int i=0; i<num_particles; i++) {
    
        old_x     = particles[i].x;
        old_y     = particles[i].y;
        old_theta = particles[i].theta;
        
        if (yaw_rate_not_zero) {
        
            // Bicycle motion model from lecture
            v_div_yawrate = velocity/yaw_rate;
            new_theta = old_theta + yaw_rate*delta_t;
            new_x = old_x + (v_div_yawrate * (sin(new_theta)-sin(old_theta)));
            new_y = old_y + (v_div_yawrate * (cos(old_theta)-cos(new_theta)));
            
        } else {
            
            // Bicycle motion model from lecture
            new_theta = old_theta;
            new_x = old_x + (velocity*delta_t)*cos(old_theta);
            new_y = old_y + (velocity*delta_t)*sin(old_theta);
            
        } // if
        
        // Noise component
        noise_x     = dist_x(gen);
        noise_y     = dist_y(gen);
        noise_theta = dist_theta(gen);
        
        // Update prediction
        particles[i].x      = new_x + dist_x(gen);
        particles[i].y      = new_y +  dist_y(gen);
        particles[i].theta  = new_theta + dist_theta(gen);
        
        if (PDEBUG) cout << "PF:Pred #" << i+1 << " " << particles[i].id << " " << particles[i].x << " "
                                                      << particles[i].y  << " " << particles[i].theta << endl;
    } // for num_particles
} // predicition


//-----
// Step #3 - updateWeights
//
// Updates the weights of each particle using a mult-variate Gaussian distribution. The landmark observations are given
// in the VEHICLE'S coordinate system and each particle (vehicle) is located according to the MAP's coordinate system.
// Note: When I display ID index it's just the number & when I use it in index SUBTRACT 1 for arrays that index from zero!
//
// Input:
//  @param sensor_range   x Initial x position [m] (simulated estimate from GPS)
//  @param std_landmark[] std[] Array of dimension 3 [standard deviation of x [m], std of y [m], std of yaw [rad]]
//  @param observations   (of landmarks) (Vehicle space)
//  @param map_landmarks    map_landmarks               (Map Space)
//
// Using:
// Particles    p                           (Map Space)
//
// Updates???
// LandmarkObs observations_onmap           (Map Space)
//-----
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
    
    LandmarkObs lm,lm_onmap;
    Particle    p;
    vector <LandmarkObs> observations_onmap;
    vector <LandmarkObs> observations_onmap_inrange;
    vector <LandmarkObs> map_landmarks_in_range;
    
    int index_into_landmark;
    double lmx, lmy, distance;
   
    double costheta,sintheta;
    double mu_x, mu_y, x_diff, y_diff, x_diff_sq, y_diff_sq;
    double sigma_x, sigma_y, sigmax_sq, sigmay_sq;
    double cterm, prob, total_prob;
    

    // Pre-compute fixed prob terms
    sigma_x = std_landmark[0];
    sigma_y = std_landmark[1];
    sigmax_sq = sigma_x*sigma_x;
    sigmay_sq = sigma_y*sigma_y;
    cterm = 1.0/(2.0*M_PI*sigma_x*sigma_y);
    
    
    // Step 1: Loop over all particles
    for (int i=0; i<particles.size(); i++) {

        p.id     = particles[i].id;
        p.x      = particles[i].x;
        p.y      = particles[i].y;
        p.theta  = particles[i].theta;
        p.weight = particles[i].weight;
        
        // Step 1A: Convert vehicle obs to map coord (from particle perspective which is in Map coord)
        observations_onmap.clear();
        observations_onmap_inrange.clear();
        for (int j=0; j<observations.size(); j++) {
        
            costheta = cos(p.theta);
            sintheta = sin(p.theta);
            
            lm.x = observations[j].x;
            lm.y = observations[j].y;
            
            lm_onmap.id = observations[j].id;
            lm_onmap.x  = p.x + (lm.x*costheta - lm.y*sintheta);
            lm_onmap.y  = p.y + (lm.x*sintheta + lm.y*costheta);
            observations_onmap.push_back(lm_onmap);
        }
        
        if (PDEBUG) {
            cout << "PF:UW # of oservs, observs_inrange, map landmarks =" << observations.size() << " " << observations_onmap.size() << " " << map_landmarks.landmark_list.size() <<  endl;
        }

        // Step 1B: Determine which map landmarks are in range of particle/vehicle sensor
        for (int j=0; j<map_landmarks.landmark_list.size(); j++) {
        
            lmx = (double)map_landmarks.landmark_list[j].x_f;
            lmy = (double)map_landmarks.landmark_list[j].y_f;
            distance = dist(p.x, p.y, lmx, lmy);

            if (distance <= sensor_range) {
                lm.id = map_landmarks.landmark_list[j].id_i;
                lm.x = map_landmarks.landmark_list[j].x_f;
                lm.y = map_landmarks.landmark_list[j].y_f;
                map_landmarks_in_range.push_back(lm);
            }
        
        }
        if (PDEBUG) cout << "PF:UW # of map landmarks in range=" << map_landmarks_in_range.size() << endl;
        
        
        // Step #2: Associate in range map landmarks to vehicle observations
        dataAssociation(map_landmarks_in_range, observations_onmap);
 
        // Check that association worked
        if (PDEBUG) {
            for (int k=0; k<observations_onmap.size(); k++) {
                int index_into_landmark = observations_onmap[k].id;
                cout << "PF:UW A=" << observations_onmap[k].x << " " << observations_onmap[k].y << " "
                     << index_into_landmark << " " << map_landmarks.landmark_list[index_into_landmark-1].x_f << " "
                                                   << map_landmarks.landmark_list[index_into_landmark-1].y_f << endl;
            }
        }
    
        
        // Step #3 - Calculate total probability and update weights
        total_prob = 1.0;
        for (int j=0; j<observations_onmap.size(); j++) {
            
            // Find map landmark associated with current observation
            index_into_landmark = observations_onmap[j].id;
          
            mu_x = map_landmarks.landmark_list[index_into_landmark-1].x_f;
            mu_y = map_landmarks.landmark_list[index_into_landmark-1].y_f;
           
            x_diff = observations_onmap[j].x - mu_x;
            y_diff = observations_onmap[j].y - mu_y;
            x_diff_sq = x_diff * x_diff;
            y_diff_sq = y_diff * y_diff;

            prob = cterm * exp(-((x_diff_sq/(2.0*sigmax_sq)) + (y_diff_sq/(2.0*sigmay_sq))));
            total_prob = total_prob * prob;
        }
        if (PDEBUG) cout << "PF:UW: Total_prob=" << total_prob << endl;

        // Update particle weight
        particles[i].weight = total_prob;
        weights[i] = total_prob;

    } // for particles
    
} // updateWeights


//-----
// Helper Method - dataAssociation (called by updateWeights) 
//
// Find the Map landmark (predicted measurement) that is closest to each observed (vehicle) landmark measurement and
// assign the observed measurement to this particular observed landmark measurement's ID. Right now, this is using
// nearest neighbor based on distance as the metric.
// Note: when this routine is called, observed landmarks are ASSUMED TO BE IN MAP coordinates!
//
// Using:
// LandmarkObs  predicted                   (Map Space - NOTE: This is the in-range map landmarks (not all of them)
// LandmarkObs  observations (of landmarks) (Map Space - Passed in already converted to Map Space)
//
// Output:
// LandmarkObs  observations (Passed by reference so updated with ID of closest map landdmark)
//-----
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
    
    LandmarkObs   o_lm = LandmarkObs();
    LandmarkObs map_lm = LandmarkObs();
    int         min_distance_lm_id;
    double      min_distance;


    // For each map landmark ("predicted") find the observation wich is closest (nearest neighbor)
    for (int i=0; i<observations.size(); i++) {
   
        // vehicle observed map landmark (in map coordinates)
        o_lm.x = observations[i].x;
        o_lm.y = observations[i].y;

        // Find smallest distance from map landmark to observed landmarks. Put in observation ID!
        min_distance = __DBL_MAX__;  // init each iter
        min_distance_lm_id = 0;
        for (int j=0; j<predicted.size(); j++) {
            map_lm.x = predicted[j].x;
            map_lm.y = predicted[j].y;
            double distance = dist(o_lm.x, o_lm.y, map_lm.x, map_lm.y);
            if (distance < min_distance) {
                min_distance_lm_id = predicted[j].id; // This is Map ID that has min distance
                min_distance = distance;
            }
            observations[i].id = min_distance_lm_id;  // Update observation with closest predicition ID
        } // for map landmarks in range
    
    } // for observations

} // dataAssociation


//-----
// Step #4 - resample
//
// Resample particles and with replacement with probability proportional to their weight.
// Using Class variables:
//   weights    <vector>
//   particles  <vector>
//-----
void ParticleFilter::resample() {
  
    // Use discrete_distribution (integer) using weights array which will then
    // assign index based on probability (tricky!)
    default_random_engine gen;
    discrete_distribution<int> distribution(weights.begin(), weights.end());  // Need begin/end of weights list
    
    vector<Particle> reorder_particles;                                       // temp storage for re-ordering particles
    
    // Re-order based on probability
    for(int i = 0; i < num_particles; i++){
        reorder_particles.push_back(particles[distribution(gen)]);
    }

    // Update
    particles = reorder_particles;
}


//---------------
// Helper Methods for Particle Filter Class
// provided by Udacity
//---------------
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

//
// Programming notes:
//
// For C++ V11 and higher, could use nice object iterator, such as:

// For particles
//for (auto const &p: particles) {
// Clear landmarks associated to observations
//for (auto & obs: observations) {
//    obs.id = 0;
//}
//
//for (auto const &map_landmark: map_landmarks.landmark_list) {
//    double distance = dist(p.x, p.y, map_landmark.x_f, map.landmark.y_f);
//    if (distance > sensor_range) {
