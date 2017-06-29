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
bool PDEBUG = false;

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

// 3 parts to predicition
//P-1 Sample
//P-2 Calculate weight
//P-3 Update state vector



// Section #1
void ParticleFilter::init(double gps_x, double gps_y, double init_heading, double gps_std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // NOTE: these are in map coordinates and are gaus around initial GPS
    //
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

    Particle a_particle;      // single particle object
    default_random_engine gen;  // rv sequence object
    
    
    // Set # of particles
    num_particles = 999;
    if (PDEBUG) cout << "PF-Init:Number of particles being generated=" << num_particles << endl;
    
    // Normal (Gaussian) distributions with means position and yaw w/ st deviations's of std[].
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
        
        particles.push_back(a_particle);  // Add to array of particles.
        weights.push_back(1.0F);          // List of NORMALIZED particle weights <TODO> Figure this out
        
        //
        if (PDEBUG) cout << "PF-Init:Particle #" << i+1 << " " << a_particle.id << " " << a_particle.x << " "
                                                               << a_particle.y  << " " << a_particle.theta << " " <<  endl;
    }
    
    is_initialized = true;
    
} // init


// Section #2
// Prediction step by adding measurements to each particle and add random Gaussian noise.
// std is deviations of position and heading angle measurements


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	
    Particle one_particle;      // single particle object
    default_random_engine gen;  // rv sequence object
    bool yaw_rate_not_zero;
    double old_x, new_x, old_y, new_y, old_theta, new_theta;
    double noise_x, noise_y, noise_theta;
    double v_div_yawrate;
    
 
    // Create gaussian distributions generators with zero means and measurement stds so can add to prediction
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
        
            // Equations from lecture
            v_div_yawrate = velocity/yaw_rate;
            new_theta = old_theta + yaw_rate*delta_t;
            new_x = old_x + (v_div_yawrate * (sin(new_theta)-sin(old_theta)));
            new_y = old_y + (v_div_yawrate * (cos(old_theta)-cos(new_theta)));
            
        } else {
            
            // Equations from lecture
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
    
} // predicition method


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
    
    // Using:
    // Particles    p                           (Map Space)
    // LandmarkObs  observations (of landmarks) (Vehicle space)
    // Map          map_landmarks               (Map Space)
    //
    // Create:
    // LandmarkObs observations_onmap           (Map Space)

    //double obsx_map, obsy_map;
    
    vector <LandmarkObs> map_observ_inrange;
    
    // For each particle
    //for (auto & p: particles) {
        // Clear landmarks associated to observations
        //for (auto & obs: observations) {
        //    obs.id = 0;
        //}

    // For each particle
    //cout << "updateWeights: # of particles=" << particles.size() << endl;
    for (int i=0; i<particles.size(); i++) {
    
        /*
        Particle a_particle;
        int    p_id     = a_particle.id;
        double p_x      = a_particle.x;
        double p_y      = a_particle.y;
        double p_yaw    = a_particle.theta;
        double p_weight = a_particle.weight;
        */
         
         
        Particle p;
        p.id     = particles[i].id;
        p.x      = particles[i].x;
        p.y      = particles[i].y;
        p.theta  = particles[i].theta;
        p.weight = particles[i].weight;
        
        
        // Step 1A: Convert vehicle obs to map coord (from particle perspective which is in Map coord)
        vector <LandmarkObs> observations_onmap;
        vector <LandmarkObs> observations_onmap_inrange;
        //observations_onmap = LandmarkObs();  // <TODO> <Need to INIT THIS!!>
        LandmarkObs lm,lm_onmap;
       

        if (PDEBUG) cout << "PF:UW # of map observations=" << observations.size() << endl;
        for (int j=0; j<observations.size(); j++) {
        
            double costheta = cos(p.theta);
            double sintheta = sin(p.theta);
            
            lm.x = observations[j].x;
            lm.y = observations[j].y;
            
            // <TODO> !!! Might need push here
            //observations_onmap[j].x = p.x + (lm.x*costheta - lm.y*sintheta);
            //observations_onmap[j].y = p.y + (lm.x*sintheta + lm.y*costheta);
            
            lm_onmap.id = observations[j].id;
            lm_onmap.x  = p.x + (lm.x*costheta - lm.y*sintheta);
            lm_onmap.y  = p.y + (lm.x*sintheta + lm.y*costheta);
            observations_onmap.push_back(lm_onmap);
            //cout << "T: P#" << i+1 << " OBS#" << j+1 << " 0=" << lm.x << " " << lm.y <<
            //                                            " 1=" << lm_onmap.x << " " << lm_onmap.y << " a=" << p.theta << endl;
        }
        
        
        // Step 1B: Determine which map landmarks are in range of particle/vehicle sensor
        std::vector<LandmarkObs> map_landmarks_in_range;
        LandmarkObs cur_lm;
        if (PDEBUG) cout << "PF:UW # of map landmarks=" << map_landmarks.landmark_list.size() << endl;;
        for (int j=0; j<map_landmarks.landmark_list.size(); j++) {
        
            double lmx = (double) map_landmarks.landmark_list[j].x_f;
            float lmy = map_landmarks.landmark_list[j].y_f;
            double distance = dist(p.x, p.y, lmx, lmy);

            if (distance <= sensor_range) {
                cur_lm.id = map_landmarks.landmark_list[j].id_i;
                cur_lm.x = map_landmarks.landmark_list[j].x_f;
                cur_lm.y = map_landmarks.landmark_list[j].y_f;

                map_landmarks_in_range.push_back(cur_lm);
            }
        
        }
        if (PDEBUG) cout << "PF:UW # of map landmarks in range=" << map_landmarks_in_range.size() << endl;
        
        // Change to this later to clean up
        //for (auto const &map_landmark: map_landmarks.landmark_list) {
        //    double distance = dist(p.x, p.y, map_landmark.x_f, map.landmark.y_f);
        //    if (distance > sensor_range) {
        
        
        
        
        //Keep only in-range observations
        /*
        for (int j=0; j<observations.size(); j++) {
            
            double p_to_obs_dist = dist(p.x, p.y, observations_onmap[j].x, observations_onmap[j].y);
            if (p_to_obs_dist < sensor_range) {
                observations_onmap_inrange.push_back(observations_onmap[j]);
                cout << "Keep: P#" << i+1 << " OBS#" << j+1 << " dist=" << p_to_obs_dist << endl;
            } else {
                cout << "DEL : P#" << i+1 << " OBS#" << j+1 << " dist=" << p_to_obs_dist << endl;
            }
        }
        cout << "After in range check: # of observations in kept=" << observations_onmap_inrange.size() << endl;;
         */
         
        // Step #2: Associate in range observations with map landmarks
        // Attempt to associate a landmark to each observation
        dataAssociation(map_landmarks_in_range, observations_onmap);
 
        // Check that we have right stuff
        if (PDEBUG) cout << "For obser #1=" << observations_onmap[0].x << " " << observations_onmap[0].y << " " << observations_onmap[0].id << endl;
        int index_into_landmark = observations_onmap[0].id;
        if (PDEBUG) cout << "Map Index for Obser=" <<  index_into_landmark-1 << endl;
        if (PDEBUG) cout << "Assoc Map=" << map_landmarks.landmark_list[index_into_landmark-1].x_f << " "
        << map_landmarks.landmark_list[index_into_landmark-1].y_f << endl;
        
        // Need to print out whole array
        // cout << "Map X= ";
        //for(int i=0 ; i<map_landmarks.landmark_list.size(); i++)
        //{
        //    cout << i << " " << map_landmarks.landmark_list[i].x_f << " ";
        //}
        //cout << endl;
      

        //cout << "Map Y= ";
        //for(int i=0 ; i<map_landmarks.landmark_list.size(); i++)
        //{
        //    cout << i << " " << map_landmarks.landmark_list[i].y_f << " ";
        //}
        //cout << endl;
        
        
        // Step #3 - Use the multi-variate Guass prob function to assign weight to each particle
        //    based on ....
        
        
        // In lesson 13, I do however see the following formula
        //P(x,y)= (1/2pi*sigmax*sigmay) * exp(- ( ((x - mux)^2)/(2*sigmax^2) + ((y - muy)^2)/(2*sigmay^2) ))

        
        
        double sigma_x = std_landmark[0];
        double sigma_y = std_landmark[1];
        
        double sigmax_sq = sigma_x*sigma_x;
        double sigmay_sq = sigma_y*sigma_y;
        
        double cterm = 1.0/(2.0*M_PI*sigma_x*sigma_y);
        //cout << "cterm=" << cterm << endl;

        
        double total_prob = 1.0;
        for (int j=0; j<observations_onmap.size(); j++) {
            
            // Find map landmark associated with current observation
            int index_into_landmark = observations_onmap[j].id;
            //cout << "At observ #=" << j+1 << " the id index into map=" << index_into_landmark << endl;
            
            double mu_x = map_landmarks.landmark_list[index_into_landmark-1].x_f;
            double mu_y = map_landmarks.landmark_list[index_into_landmark-1].y_f;
            //cout << "mu_x & y=" << mu_x  << " " << mu_y  << endl;
            
            double x_diff = observations_onmap[j].x - mu_x;
            double y_diff = observations_onmap[j].y - mu_y;
            //cout << "x & y diff=" << x_diff  << " " << y_diff << endl;
            
            double x_diff_sq = x_diff * x_diff;
            double y_diff_sq = y_diff * y_diff;

            double prob = cterm * exp(-((x_diff_sq/(2.0*sigmax_sq)) + (y_diff_sq/(2.0*sigmay_sq))));
            //cout << "prob=" << prob << endl;
            
            total_prob = total_prob * prob;
            //cout << "total_prob=" << total_prob << endl;
            
        }
        if (PDEBUG) cout << "PF:UW: Total_prob=" << total_prob << endl;

        
        // Update particle weight
        particles[i].weight = total_prob;
        weights[i] = particles[i].weight;

        
        /*
        // Create a list of in-range observations in map coordiantes
        // #1A - transform vehicle landmark observations to map coordinates
        cout << "updateWeights: # of observations=" << observations.size() << endl;;
        for (int j=0; j<observations.size(); j++) {
            
            LandmarkObs lm;
            lm.x = observations[j].x;
            lm.y = observations[j].y;

            //double dist = dist(p.x, p.y, landmark_s.x_f, landmark_s.y_f);
            if (PDEBUG) cout << "LM #" << observations[i].id << " x=" << lm.x << " y="
                << lm.y << endl;
            double particleto_observedlm_dist = dist(p.x, p.y, lm.x, lm.y);
            cout << "updateWeights: obs #" << j+1 << "  dist=" << particleto_observedlm_dist << endl;;
            if (particleto_observedlm_dist < sensor_range) {
                
                // Convert coords to Map
                double cost = cos(p.theta);
                double sint = sin(p.theta);
                LandmarkObs map_obs = LandmarkObs();
                map_obs.x = p.x + (lm.x * cost - lm.y * sint);
                map_obs.y = p.y + (lm.x * sint + lm.y * cost);
            
                // Push to observed in range in map coord
                map_observ_inrange.push_back(map_obs);
                cout << "updateWeights: obs in range #" << j+1 << "  x=" << map_obs.x
                        << "  y=" << map_obs.y <<endl;;
            }
          */
        // Call data association

        
        //obs_x_on_map = observations[i].x +
    
    //}
    
    /*
    // Convert observations from VEHICULE to MAP's coordinate
    std::vector<LandmarkObs> map_observations;
    for (auto const& obs: observations) {
        double cos_t = cos(p.theta);
        double sin_t = sin(p.theta);
        LandmarkObs map_obs = LandmarkObs();
        map_obs.x = p.x + (obs.x * cos_t - obs.y * sin_t);
        map_obs.y = p.y + (obs.x * sin_t + obs.y * cos_t);
        
        map_observations.push_back(map_obs);
    }
    
    // Optimization - Filter out landmarks that are out of sensor range for current particle
    std::vector<LandmarkObs> map_landmarks_in_range;
    for (auto const& landmark_s: map_landmarks.landmark_list) {
        double distance = dist(p.x, p.y, landmark_s.x_f, landmark_s.y_f);
        if (distance > sensor_range) {
            // Filter out
            continue;
        }
        
        LandmarkObs landmark = LandmarkObs();
        landmark.id = landmark_s.id_i;
        landmark.x = landmark_s.x_f;
        landmark.y = landmark_s.y_f;
        map_landmarks_in_range.push_back(landmark);
    }
    */
    
    // #1B - Add only landmarks within sensor range of
    
    
    
    // #2 - For each observed landmark , predict which landmark it goes with (within sensor range) by nearest neighbor
    
  
    } // for particles
    
    /*
    // Final Step - Normalize particle weights
    double weight_sum = 0.0F;
    for (int i=0; i<particles.size(); i++) {
        weight_sum += particles[i].weight;
    }
    for (int i=0; i<particles.size(); i++) {
        particles[i].weight =   particles[i].weight/weight_sum;
    }
    */
    
} // updateWeights


// Called by updateWeights
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.
    
// Obs are in map space at this point

    // Using:
    // LandmarkObs  predicted                   (Map Space) - This is the in-range map landmarks
    // LandmarkObs  observations (of landmarks) (Passed in converted to Map Space)
    //
    // Create:
    // LandmarkObs observations (of lan          (Map Space)
    //
    // Output:
    // LandmarkObs predicited (Map Space)

    
    
    
    // For each Map Landmark ("predicted") find the observation with is closest (nearest neighbor)
    for (int i=0; i<observations.size(); i++) {
   
        int min_distance_lm_id;
        double min_distance;
        // Access map landmark
        LandmarkObs o_lm = LandmarkObs();
       
        
        //map_lm.x = predicted[i].x;
        //map_lm.y = predicted[i].y;
        o_lm.x = observations[i].x;
        o_lm.y = observations[i].y;

        // Find smallest distance from map landmark to observed landmarks
        // Observation needs the prediction id!
        min_distance = __DBL_MAX__;
        min_distance_lm_id = 0;
        
        LandmarkObs min_map_lm = LandmarkObs();

        
        for (int j=0; j<predicted.size(); j++) {
            LandmarkObs map_lm = LandmarkObs();
                      //o_lm.id = observations[j].id;
            map_lm.x = predicted[j].x;
            map_lm.y = predicted[j].y;
            double distance = dist(o_lm.x, o_lm.y, map_lm.x, map_lm.y);
            //cout << "OBS#=" << i+1 << " MAP#=" << j+1 << " dist=" << distance << endl;
            if (distance < min_distance) {
                min_distance_lm_id = predicted[j].id; // This is Map ID that has min distance
                min_distance = distance;
                min_map_lm.x = map_lm.x;
                min_map_lm.y = map_lm.y;
                //cout << "New min dist=" << min_distance << " " << min_distance_lm_id << endl;
                ////cout << "Map coords=" << predicted[min_distance_lm_id].x << " " << predicted[min_distance_lm_id].y << endl;
                //cout << "Map coords=" << min_map_lm.x << " " << min_map_lm.y << endl;
            }
            
        //predicted[i].id = min_distance_lm_id; // Update predicition ID with the closest observation id
        observations[i].id = min_distance_lm_id;
        }
        
        //cout << "OBS#=" << i+1 << " gets map lm id#=" << min_distance_lm_id << endl;
        //cout << "OBS coords=" << o_lm.x << " " << o_lm.y << endl;
        //cout << "Nearest Map coords=" <<  min_map_lm.x  << " " <<  min_map_lm.y  << endl;
            
            
        
        // For min distance, assign to precited landmarks this observation id
    }
} // dataAssociation



// Section #4
void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    default_random_engine gen;
    discrete_distribution<int> distribution(weights.begin(), weights.end());
    
    vector<Particle> resample_particles;
    for(int i = 0; i < num_particles; i++){
        resample_particles.push_back(particles[distribution(gen)]);
        //cout << resample_particles[i].weight << " ";
    }
    //cout << endl;
    
    particles = resample_particles;
}









//---------------
// Helper Methods for Particle Filter class programmed by Udacity
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
