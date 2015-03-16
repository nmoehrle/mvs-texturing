#include "texturing.h"

void
run_mrf_optimization(MRF * mrf, const Arguments &conf) {
    util::WallTimer timer;

    std::vector<MRF::ENERGY_TYPE> energies;

    std::cout << "\tIteration\tEnergy\tRuntime" << std::endl;
    MRF::ENERGY_TYPE const zero = MRF::ENERGY_TYPE(0);
    MRF::ENERGY_TYPE last_energy = zero;
    MRF::ENERGY_TYPE energy = mrf->compute_energy();
    MRF::ENERGY_TYPE diff = last_energy - energy;

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

    if (conf.write_mrf_energies)
        write_vector_to_csv(conf.out_prefix + "_mrf_energies.csv", energies, "Energy");

    if (diff == zero)
        std::cout << "\t" << "Converged" << std::endl;
    if (diff < zero)
        std::cout << "\t" << "Increase of energy detected! Aborting..." << std::endl;
    std::cout << "\tTook: " << timer.get_elapsed_sec() << "s" << std::endl;
}
