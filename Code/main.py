# Utilities
from utilities import *


if __name__ == "__main__":

    # Input files
    file_in_NBC_correctors = "NBC_correctors.txt"
    file_in_OBC_correctors = "OBC_correctors.txt"
    file_in_NBCCYC_correctors = "NBCCYC_correctors.txt"
    file_in_OBCCYC_correctors = "OBCCYC_correctors.txt"

    target_in_min_entropy_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    # Output files
    file_out_NBC_appropriate_correctors = 'NBC_Appropriate_correctors.txt'
    file_out_OBC_appropriate_correctors = 'OBC_Appropriate_correctors.txt'
    file_out_NBCCYC_appropriate_correctors = 'NBCCYC_Appropriate_correctors.txt'
    file_out_OBCCYC_appropriate_correctors = 'OBCCYC_Appropriate_correctors.txt'
    file_out_hmin_diff = 'Hmin_in_differences.txt'
    file_out_write_OBC_optimal = 'Optimal_correctors_OBC.txt'
    file_out_write_NBC_optimal = 'Optimal_correctors_NBC.txt'
    file_out_write_OBCCYC_optimal = 'Optimal_correctors_OBCCYC.txt'
    file_out_write_NBCCYC_optimal = 'Optimal_correctors_NBCCYC.txt'

    file_out_write_targeted_Pareto_NBC = 'Optimal_correctors_NBC_9_targeted_H_in.txt'
    file_out_write_targeted_Pareto_NBCCYC = 'Optimal_correctors_NBCCYC_9_targeted_H_in.txt'
    file_out_write_targeted_Pareto_OBC = 'Optimal_correctors_OBC_9_targeted_H_in.txt'
    file_out_write_targeted_Pareto_OBCCYC = 'Optimal_correctors_OBCCYC_9_targeted_H_in.txt'

    # Output lists
    NBC_entropy_bound_code_list = list()
    OBC_entropy_bound_code_list = list()
    NBCCYC_entropy_bound_code_list = list()
    OBCCYC_entropy_bound_code_list = list()
    code_weight_distribution_list = list()
    differences_OBC_and_NBC_bound_list = list()
    optimal_correctors_for_targeted_h_in_list_from_OBC = list()
    optimal_correctors_for_targeted_h_in_list_from_OBCCYC = list()
    optimal_correctors_for_targeted_h_in_list_from_NBC = list()
    optimal_correctors_for_targeted_h_in_list_from_NBCCYC = list()

    # Load the weight distributions
    print("Load the weight distributions...")
    load_weight_distributions (file_in_NBC_correctors, code_weight_distribution_list)
    # Find appropriate correctors in OBC
    print("Find appropriate correctors in OBC...")
    find_appropriate_correctors_in_OBC (OBC_entropy_bound_code_list, file_in_OBC_correctors, file_out_OBC_appropriate_correctors, 'OBC')
    # Find appropriate correctors in OBCCYC
    print("Find appropriate correctors in OBCCYC...")
    find_appropriate_correctors_in_OBC (OBCCYC_entropy_bound_code_list, file_in_OBCCYC_correctors, file_out_OBCCYC_appropriate_correctors, 'OBCCYC')
    # Find appropriate correctors in NBC
    print("Find appropriate correctors in NBC...")
    find_appropriate_correctors_in_NBC (NBC_entropy_bound_code_list, file_in_NBC_correctors, file_out_NBC_appropriate_correctors, 'NBC')
    # Find appropriate correctors in NBCCYC
    print("Find appropriate correctors in NBCCYC...")
    find_appropriate_correctors_in_NBC (NBCCYC_entropy_bound_code_list, file_in_NBCCYC_correctors, file_out_NBCCYC_appropriate_correctors, 'NBCCYC')

    # Determine the differences between the NBC and OBC bounds
    print("Determine the differences between the NBC and OBC bounds...")
    sorted_list_codes_by_OBC_NBC_difference, sorted_list_codes_by_OBC_NBC_relative_difference = determine_bound_differences (file_out_hmin_diff, NBC_entropy_bound_code_list, OBC_entropy_bound_code_list, differences_OBC_and_NBC_bound_list)


    # Find optimal correctors from the appropriate correctors in OBC
    print("Find optimal correctors from the appropriate correctors in OBC...")
    Optimal_OBC_correctors_code_rate, Optimal_OBC_correctors_in_min_entropy, Optimal_OBC_correctors = determine_Optimal_correctors (OBC_entropy_bound_code_list, file_out_write_OBC_optimal, 'OBC')
    # Find optimal correctors from the appropriate correctors in NBC
    print("Find optimal correctors from the appropriate correctors in NBC...")
    Optimal_NBC_correctors_code_rate, Optimal_NBC_correctors_in_min_entropy, Optimal_NBC_correctors = determine_Optimal_correctors (NBC_entropy_bound_code_list, file_out_write_NBC_optimal, 'NBC')
    # Find optimal correctors from the appropriate correctors in OBCCYC
    print("Find optimal correctors from the appropriate correctors in OBCCYC...")
    Optimal_OBCCYC_correctors_code_rate, Optimal_OBCCYC_correctors_in_min_entropy, Optimal_OBCCYC_correctors = determine_Optimal_correctors (OBCCYC_entropy_bound_code_list, file_out_write_OBCCYC_optimal, 'OBCCYC')
    # Find optimal correctors from the appropriate correctors in NBCCYC
    print("Find optimal correctors from the appropriate correctors in NBCCYC...")
    Optimal_NBCCYC_correctors_code_rate, Optimal_NBCCYC_correctors_in_min_entropy, Optimal_NBCCYC_correctors = determine_Optimal_correctors (NBCCYC_entropy_bound_code_list, file_out_write_NBCCYC_optimal, 'NBCCYC')

    # Find optimal correctors for targeted h_in in OBC
    print("Find optimal correctors for targeted h_in in OBC...")
    determine_Optimal_correctors_for_targeted_h_in (Optimal_OBC_correctors, file_out_write_targeted_Pareto_OBC, optimal_correctors_for_targeted_h_in_list_from_OBC, 'OBC', target_in_min_entropy_list)
    # Find optimal correctors for targeted h_in in OBCCYC
    print("Find optimal correctors for targeted h_in in OBCCYC...")
    determine_Optimal_correctors_for_targeted_h_in (Optimal_OBCCYC_correctors, file_out_write_targeted_Pareto_OBCCYC, optimal_correctors_for_targeted_h_in_list_from_OBCCYC, 'OBCCYC', target_in_min_entropy_list)
    # Find optimal correctors for targeted h_in in NBC
    print("Find optimal correctors for targeted h_in in NBC...")
    determine_Optimal_correctors_for_targeted_h_in (Optimal_NBC_correctors, file_out_write_targeted_Pareto_NBC, optimal_correctors_for_targeted_h_in_list_from_NBC, 'NBC', target_in_min_entropy_list, code_weight_distribution_list)
    # Find optimal correctors for targeted h_in in NBCCYC
    print("Find optimal correctors for targeted h_in in NBCCYC...")
    determine_Optimal_correctors_for_targeted_h_in (Optimal_NBCCYC_correctors, file_out_write_targeted_Pareto_NBCCYC, optimal_correctors_for_targeted_h_in_list_from_NBCCYC, 'NBCCYC', target_in_min_entropy_list, code_weight_distribution_list)

    # Find the extraction efficiency for the entire range of h_in
    print("Find the extraction efficiency for the entire range of h_in...")
    OBC_in_min_entropies_range, Optimal_OBC_efficiencies, NBC_in_min_entropies_range, Optimal_NBC_efficiencies = find_efficiencies_for_entire_h_in_range (Optimal_OBC_correctors, Optimal_NBC_correctors, code_weight_distribution_list)
    OBCCYC_in_min_entropies_range, Optimal_OBCCYC_efficiencies, NBCCYC_in_min_entropies_range, Optimal_NBCCYC_efficiencies = find_efficiencies_for_entire_h_in_range (Optimal_OBCCYC_correctors, Optimal_NBCCYC_correctors, code_weight_distribution_list)

    # Generate plots
    print("Generate plots...")

    graph_abs_and_rel_entropy_differences_OBC_and_NBC_title = 'Abs_and_rel_entropy_differences_OBC_and_NBC.pdf'
    PLOT_abs_and_rel_entropy_differences_OBC_and_NBC (sorted_list_codes_by_OBC_NBC_difference, sorted_list_codes_by_OBC_NBC_relative_difference, code_weight_distribution_list, graph_abs_and_rel_entropy_differences_OBC_and_NBC_title)

    graph_OBC_NBC_correctors_and_efficiencies_title = 'OBC_NBC_correctors_and_efficiencies.pdf'
    PLOT_OBC_NBC_correctors_and_efficiencies (OBC_entropy_bound_code_list, NBC_entropy_bound_code_list, Optimal_OBC_correctors_in_min_entropy, Optimal_NBC_correctors_in_min_entropy, Optimal_OBC_correctors_code_rate,  Optimal_NBC_correctors_code_rate, OBC_in_min_entropies_range, Optimal_OBC_efficiencies, NBC_in_min_entropies_range, Optimal_NBC_efficiencies, graph_OBC_NBC_correctors_and_efficiencies_title)

    graph_OBCCYC_NBCCYC_correctors_title = 'OBCCYC_NBCCYC_correctors.pdf'
    PLOT_OBCCYC_NBCCYC_correctors (Optimal_OBCCYC_correctors_in_min_entropy, Optimal_NBCCYC_correctors_in_min_entropy, Optimal_OBCCYC_correctors_code_rate,  Optimal_NBCCYC_correctors_code_rate, graph_OBCCYC_NBCCYC_correctors_title)

    graph_OBCCYC_NBCCYC_efficiencies_title = 'OBCCYC_NBCCYC_efficiencies.pdf'
    PLOT_OBCCYC_NBCCYC_efficiencies (OBCCYC_in_min_entropies_range, Optimal_OBCCYC_efficiencies, NBCCYC_in_min_entropies_range, Optimal_NBCCYC_efficiencies, graph_OBCCYC_NBCCYC_efficiencies_title)
