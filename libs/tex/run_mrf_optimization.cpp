#include "texturing.h"

void
run_mrf_optimization(mrf::Graph * mrf) {
    util::WallTimer timer;

    std::vector<mrf::ENERGY_TYPE> energies;

    std::cout << "\tIteration\tEnergy\tRuntime" << std::endl;
    mrf::ENERGY_TYPE const zero = mrf::ENERGY_TYPE(0);
    mrf::ENERGY_TYPE last_energy = zero;
    mrf::ENERGY_TYPE energy = mrf->compute_energy();
    mrf::ENERGY_TYPE diff = last_energy - energy;

    unsigned int i = 0;
    while (diff != zero) {
        std::cout << "\t" << i << "\t" << energy << "\t" << timer.get_elapsed_sec() << std::endl;
        energies.push_back(energy);
        last_energy = energy;
        i++;
        energy = mrf->optimize(1);
        diff = last_energy - energy;
        if (diff <= 0)
            break;
    }
    std::cout << "\t" << i << "\t" << energy << std::endl;

    //if (conf.write_mrf_energies)
    //    write_vector_to_csv(conf.out_prefix + "_mrf_energies.csv", energies, "Energy");

    if (diff == zero)
        std::cout << "\t" << "Converged" << std::endl;
    if (diff < zero)
        std::cout << "\t" << "Increase of energy detected! Aborting..." << std::endl;
    std::cout << "\tTook: " << timer.get_elapsed_sec() << "s" << std::endl;
}
