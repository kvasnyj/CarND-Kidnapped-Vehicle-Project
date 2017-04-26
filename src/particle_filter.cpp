/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"
#include "helper_functions.h"

std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).


    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_psi(theta, std[2]);

    num_particles = 1000;
    for (int i = 0; i < num_particles; i++) {
        Particle p =
                {
                        i, // id
                        dist_x(gen), // x
                        dist_y(gen), // y
                        dist_psi(gen), // theta
                        1 //weight
                };

        particles.push_back(p);
    }

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    for (int i = 0; i < num_particles; i++) {
        Particle p = particles[i];

        if (yaw_rate == 0) // going straight
        {
            p.x += velocity * delta_t * sin(p.theta);
            p.y += velocity * delta_t * cos(p.theta);
        } else {
            p.x += velocity / yaw_rate * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
            p.y += velocity / yaw_rate * (cos(p.theta + yaw_rate * delta_t) - cos(p.theta));
            p.theta += yaw_rate * delta_t;
        }


    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

    for (int i = 0; i < observations.size(); i++) {
        double x1 = observations[i].x;
        double y1 = observations[i].y;

        double min_d = -1;
        for (int j = 0; j < predicted.size(); j++) {
            double x2 = predicted[j].x;
            double y2 = predicted[j].y;

            double d = dist(x1, y1, x2, y2);
            if (min_d > d || min_d < 0) {
                observations[i].id = j;
                min_d = d;
            }
        }
    }

}

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
    //   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
    //   for the fact that the map's y-axis actually points downwards.)
    //   http://planning.cs.uiuc.edu/node99.html


    //dataAssociation(predicted_landmarks, observations);
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

void ParticleFilter::write(std::string filename) {
    // You don't need to modify this file.
    std::ofstream dataFile;
    dataFile.open(filename, std::ios::app);
    for (int i = 0; i < num_particles; ++i) {
        dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
    }
    dataFile.close();
}
