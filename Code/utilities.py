from parameters import *

# %%
def set_size(width, golden_ratio_boolean=True, new_golden_ratio=0.5, fraction=1, subplots=(1, 1), width_in_pt = True):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float
            Document textwidth or columnwidth in pts
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width * fraction

    # Convert from pt to inches
    if (width_in_pt):
        inches_per_pt = 1 / 72.27
    else:
        inches_per_pt = 1

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    if (golden_ratio_boolean):
        golden_ratio = (5**.5 - 1) / 2
    else:
        golden_ratio = new_golden_ratio

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim

# %%
"""Determine codes on the Optimal Frontier - 2D"""
def Optimal_frontier_2d_function(Rates, Entropies, ns, ks, ds, effs, maxX = True, maxY = True):
    myList = sorted([[Rates[i], Entropies[i], ns[i], ks[i], ds[i], effs[i]] for i in range(len(Rates))], reverse=maxX)
    p_front = [myList[0]]    
    for pair in myList[1:]:
        if maxY: 
            if pair[1] > p_front[-1][1]:
                p_front.append(pair)
        else:
            if pair[1] < p_front[-1][1]:
                p_front.append(pair)
    p_front_Rates = [pair[0] for pair in p_front]
    p_front_Entropies = [pair[1] for pair in p_front]
    p_front_n = [pair[2] for pair in p_front]
    p_front_k = [pair[3] for pair in p_front]
    p_front_d = [pair[4] for pair in p_front]
    p_front_eff = [pair[5] for pair in p_front]
    return p_front_Rates, p_front_Entropies, p_front_n, p_front_k, p_front_d, p_front_eff, p_front

# %%
'''New entropy bound based on the dual weight distributions'''
def NEW_H_min_bit_bound(n, k, Hin, weight_distribution):
    """
    Calculate H_min_out using the Decimal data type.
    """
    if Hin == Decimal(0):
        return Decimal(0)
    elif Hin == Decimal(1):
        return Decimal(k)
    else:
        p_max = Decimal(0)
        two_pow_neg_k = Decimal(2) ** Decimal(-k)
        factor = (Decimal(2) ** (Decimal(1) - Hin)) - Decimal(1)

        for elem in weight_distribution:
            weight, prob = elem
            p_max += two_pow_neg_k * prob * (factor ** Decimal(weight))

        if p_max > Decimal(1) or p_max < Decimal(0):
            raise ValueError("p_max is not in (0, 1) range!")
            
        return -(p_max.log10() / Decimal(2).log10())



# %%
def bisection_NEW_Hin_higher_than_target(n_code, k_code, target_H_min_tot, code_weight_distribution, a, b, tol, max_iter):
    """
    Uses the bisection method with high precision to find the value of Hin that results in an H_min_tot slightly higher than the target.

    :param n_code: Length of the code.
    :param k_code: A parameter related to the code.
    :param target_H_min_tot: The target H_min_tot value.
    :param code_weight_distribution: List of tuples with weight and its frequency.
    :param a: Lower bound for bisection.
    :param b: Upper bound for bisection.
    :param tol: Tolerance level.
    :param max_iter: Maximum number of iterations.

    :return: A tuple containing the value of Hin that gives an H_min_tot slightly higher than the target (up to the required tolerance level) and the corresponding difference.
    """

    # Function to compute the difference between the target and computed H_min_tot
    def difference(Hin):
      return NEW_H_min_bit_bound(n_code, k_code, Hin, code_weight_distribution) - target_H_min_tot

    # Check if the provided interval is valid
    if difference(a) * difference(b) > 0:
        return Decimal(1), Decimal(1) # The requested input min-entropy value is higher than 0.999

    # Bisection algorithm
    for _ in range(max_iter):
        midpoint = (a + b) / Decimal(2)
        diff = difference(midpoint)

        # If difference is very close to zero or the interval is smaller than the tolerance
        if abs(diff) < tol or (b - a) / Decimal(2) < tol:
            # If the difference is positive, we have a candidate that gives a higher H_min_tot
            if diff > Decimal(0):
                diff_in = NEW_H_min_bit_bound(n_code, k_code, midpoint, code_weight_distribution) - target_H_min_tot
                return midpoint, diff_in
            else:
                # If the difference is negative, we continue searching for a higher H_min_tot
                a = midpoint

        # Adjust the interval for the next iteration
        elif (diff > Decimal(0) and difference(a) > Decimal(0)) or (diff < Decimal(0) and difference(a) < Decimal(0)):
            a = midpoint
        else:
            b = midpoint
    raise ValueError("Bisection method did not converge.")

# %%
'''Lacharme's (old) entropy bound'''
def OLD_H_min_bit_bound(k_code, Hin, d_code):
    p_base = Decimal(1) + (( Decimal(2)**(Decimal(1) - Hin) - Decimal(1))**d_code) * (Decimal(2)**k_code)
    H_min_tot = k_code - p_base.log10() / Decimal(2).log10()
    return H_min_tot

# %%
'''Inversed Lacharme's (old) entropy bound'''
def OLD_inverse_H_min_bit_bound(k_code, d_code, minimum_output_entropy_per_bit_bound):
    p_base = ( Decimal(2)**(-k_code/d_code) ) * ( ( Decimal(2)**(Decimal(1) - minimum_output_entropy_per_bit_bound) - Decimal(1))**(Decimal(1)/d_code) ) + Decimal(1) 
    Hin = Decimal(1) - (p_base.log10() / Decimal(2).log10())
    return Hin

# %%
'''
Applicable to: all codes
Bounds: Old and New bound
'''
def find_appropriate_correctors_in_OBC(OBC_entropy_bound_code_list, file_in_OBC_correctors, file_out_OBC_appropriate_correctors, set_name):
    total_num_codes_OBC = 0
    with open(file_in_OBC_correctors, 'r') as f_reader, open(file_out_OBC_appropriate_correctors, 'w') as f_OBC_writer:
        read_line = f_reader.readline()
        while read_line:
            if 'Code' in read_line:
                start_bracket_pos = read_line.find('[')
                end_bracket_pos = read_line.find(']')
                code_text = read_line[start_bracket_pos+1:end_bracket_pos]
                curr_code = [int(i) for i in code_text.split(',')]
                total_num_codes_OBC += 1
                output_entropy_per_bit_bound = curr_code[1] - 1 + minimum_output_entropy_per_bit_bound
                H_min_IN = OLD_inverse_H_min_bit_bound(Decimal(curr_code[1]), Decimal(curr_code[2]), minimum_output_entropy_per_bit_bound)
                if H_min_IN < minimum_output_entropy_per_bit_bound:
                    efficiency = OLD_H_min_bit_bound(curr_code[1], H_min_IN, curr_code[2]) / (H_min_IN * curr_code[0])
                    efficiency_str = str(efficiency.quantize(Decimal('.00000001'), rounding = ROUND_HALF_EVEN))# Regular rounding
                    OBC_entropy_bound_code_list.append((curr_code[0],curr_code[1], curr_code[2], curr_code[1]/curr_code[0], H_min_IN, efficiency))
                    H_in_req_str = str(H_min_IN.quantize(Decimal('.00000001'), rounding = ROUND_UP))# Conservative rounding
                    f_OBC_writer.write((f'n = {curr_code[0]}, k = {curr_code[1]}, d = {curr_code[2]}, CR = {curr_code[1]/curr_code[0]:.8f}, H_in_req = {H_in_req_str}, efficiency = {efficiency_str} \n'))
            read_line = f_reader.readline()
    
    print(f'({set_name}) Number of correctors in {set_name}: {total_num_codes_OBC}')
    print(f'({set_name}) Number of appropriate correctors in {set_name}: {len(OBC_entropy_bound_code_list)}')
    with open(file_out_OBC_appropriate_correctors, 'a') as f_OBC_writer:
        f_OBC_writer.write(f'\nNumber of correctors in {set_name}: {total_num_codes_OBC}')
        f_OBC_writer.write(f'\nNumber of appropriate correctors in {set_name}: {len(OBC_entropy_bound_code_list)}')

def find_appropriate_correctors_in_NBC(NBC_entropy_bound_code_list, file_in_NBC_correctors, file_out_NBC_appropriate_correctors, set_name):
    total_num_codes_NBC = 0
    with open(file_in_NBC_correctors, 'r') as f_reader, open(file_out_NBC_appropriate_correctors, 'w') as f_NBC_writer:
        read_line = f_reader.readline()
        while read_line:
            if 'Code' in read_line:
                start_bracket_pos = read_line.find('[')
                end_bracket_pos = read_line.find(']')
                code_text = read_line[start_bracket_pos+1:end_bracket_pos]
                curr_code = [int(i) for i in code_text.split(',')]
                weight_distribution_string = f_reader.readline()
                total_num_codes_NBC += 1
                code_weight_distribution_str = re.findall('\<[^>]*\>',weight_distribution_string)
                code_weight_distribution = [tuple(int(i) for i in el.strip('<>').split(',')) for el in code_weight_distribution_str]
                target_H_min_tot = curr_code[1] - 1 + minimum_output_entropy_per_bit_bound
                H_min_IN, diff_in = bisection_NEW_Hin_higher_than_target(curr_code[0], curr_code[1], Decimal(target_H_min_tot), code_weight_distribution, Decimal(0), Decimal('0.999'), tol, max_bisection_iter)
                if H_min_IN < minimum_output_entropy_per_bit_bound:
                    efficiency = NEW_H_min_bit_bound (curr_code[0], curr_code[1], H_min_IN, code_weight_distribution) / (H_min_IN * curr_code[0])
                    efficiency_str = str(efficiency.quantize(Decimal('.00000001'), rounding = ROUND_HALF_EVEN)) # Regular rounding
                    NBC_entropy_bound_code_list.append((curr_code[0],curr_code[1], curr_code[2], curr_code[1]/curr_code[0], H_min_IN, efficiency))           
                    H_in_req_str = str(H_min_IN.quantize(Decimal('.00000001'), rounding = ROUND_UP)) # Conservative rounding
                    f_NBC_writer.write((f'n = {curr_code[0]}, k = {curr_code[1]}, d = {curr_code[2]}, CR = {curr_code[1]/curr_code[0]:.8f}, H_in_req = {H_in_req_str}, efficiency = {efficiency_str} \n'))
            read_line = f_reader.readline()

    print(f'({set_name}) Number of correctors in {set_name}: {total_num_codes_NBC}')
    print(f'({set_name}) Number of appropriate correctors in {set_name}: {len(NBC_entropy_bound_code_list)}')
    with open(file_out_NBC_appropriate_correctors, 'a') as f_NBC_writer:
        f_NBC_writer.write(f'\nNumber of correctors in {set_name}: {total_num_codes_NBC}')
        f_NBC_writer.write(f'\nNumber of appropriate correctors in {set_name}: {len(NBC_entropy_bound_code_list)}')


def load_weight_distributions(file_in_NBC_correctors, code_weight_distribution_list):
    with open(file_in_NBC_correctors, 'r') as f_reader:
        read_line = f_reader.readline()
        while read_line:
            if 'Code' in read_line:
                start_bracket_pos = read_line.find('[')
                end_bracket_pos = read_line.find(']')
                code_text = read_line[start_bracket_pos+1:end_bracket_pos]
                curr_code = [int(i) for i in code_text.split(',')]
                weight_distribution_string = f_reader.readline()
                code_weight_distribution_str = re.findall('\<[^>]*\>',weight_distribution_string)
                code_weight_distribution = [tuple(int(i) for i in el.strip('<>').split(',')) for el in code_weight_distribution_str]
                code_weight_distribution_list.append((curr_code,code_weight_distribution))
            read_line = f_reader.readline()


# %%
'''
Applicable to: All codes
Bounds: Lacharme and new bound 
'''

def determine_bound_differences(file_name_hmin_diff, NBC_entropy_bound_code_list, OBC_entropy_bound_code_list, Differences_OBC_and_NBC_bound):

    with open(file_name_hmin_diff, 'w') as f_hmin_diff:
        for NBC_bound_code in  NBC_entropy_bound_code_list:
            found = 0
            for OBC_bound_code in OBC_entropy_bound_code_list:
                if NBC_bound_code[0] == OBC_bound_code[0] and NBC_bound_code[1] == OBC_bound_code[1] and NBC_bound_code[2] == OBC_bound_code[2]:
                    hmin_in_diff = OBC_bound_code[4] - NBC_bound_code[4]
                    hmin_in_diff_str = str(hmin_in_diff.quantize(Decimal('.00000001'), rounding = ROUND_HALF_EVEN)) # Regular rounding
                    relative_hmin_in_diff =  1 - NBC_bound_code[4]/OBC_bound_code[4]
                    relative_hmin_in_diff_str = str(relative_hmin_in_diff.quantize(Decimal('.00000001'), rounding = ROUND_HALF_EVEN)) # Regular rounding
                    Differences_OBC_and_NBC_bound.append((NBC_bound_code[0], NBC_bound_code[1], NBC_bound_code[2], NBC_bound_code[1]/NBC_bound_code[0], hmin_in_diff, relative_hmin_in_diff, OBC_bound_code[4], NBC_bound_code[4]))
                    f_hmin_diff.write((f'n = {NBC_bound_code[0]}, k = {NBC_bound_code[1]}, d = {NBC_bound_code[2]}, CR = {NBC_bound_code[1]/NBC_bound_code[0]:.8f}, Hmin Old - New bound  = {hmin_in_diff_str}, Relative Hmin improvement: {relative_hmin_in_diff_str}\n'))
                    break

    sorted_list_codes_by_OBC_NBC_difference = sorted(Differences_OBC_and_NBC_bound, key = lambda x: x[4], reverse = True)
    highest_ent_diff_str = str(sorted_list_codes_by_OBC_NBC_difference[0][4].quantize(Decimal('.00000001'), rounding = ROUND_HALF_EVEN)) # Regular rounding
    higher_h_ent_diff_str = str(sorted_list_codes_by_OBC_NBC_difference[0][6].quantize(Decimal('.00000001'), rounding = ROUND_HALF_EVEN)) # Regular rounding
    lower_h_ent_diff_str = str(sorted_list_codes_by_OBC_NBC_difference[0][7].quantize(Decimal('.00000001'), rounding = ROUND_HALF_EVEN)) # Regular rounding
    lowest_ent_diff_str = str(sorted_list_codes_by_OBC_NBC_difference[-1][4].quantize(Decimal('.00000001'), rounding = ROUND_HALF_EVEN)) # Regular rounding
    print(f'Code with the highest input min-entropy difference between OBC and NBC is code [{sorted_list_codes_by_OBC_NBC_difference[0][0]}, {sorted_list_codes_by_OBC_NBC_difference[0][1]}, {sorted_list_codes_by_OBC_NBC_difference[0][2]}] with difference of {highest_ent_diff_str} = {higher_h_ent_diff_str} - {lower_h_ent_diff_str}')
    print(f'Code with the lowest input min-entropy difference between OBC and NBC is code [{sorted_list_codes_by_OBC_NBC_difference[-1][0]}, {sorted_list_codes_by_OBC_NBC_difference[-1][1]}, {sorted_list_codes_by_OBC_NBC_difference[-1][2]}] with difference of {lowest_ent_diff_str}')

    sorted_list_codes_by_OBC_NBC_relative_difference = sorted(Differences_OBC_and_NBC_bound, key = lambda x: x[5], reverse = True)

    print(f'Code with the highest input min-entropy relative reduction between OBC and NBC is code [{sorted_list_codes_by_OBC_NBC_relative_difference[0][0]}, {sorted_list_codes_by_OBC_NBC_relative_difference[0][1]}, {sorted_list_codes_by_OBC_NBC_relative_difference[0][2]}] with reduction of {sorted_list_codes_by_OBC_NBC_relative_difference[0][5]*100:.2f} %')
    print(f'Code with the lowest input min-entropy relative reduction between OBC and NBC is code [{sorted_list_codes_by_OBC_NBC_relative_difference[-1][0]}, {sorted_list_codes_by_OBC_NBC_relative_difference[-1][1]}, {sorted_list_codes_by_OBC_NBC_relative_difference[-1][2]}] with reduction of {sorted_list_codes_by_OBC_NBC_relative_difference[-1][5]*100:.2f} %')

    with open(file_name_hmin_diff, 'a') as f_hmin_diff:
        f_hmin_diff.write(f'\nCode with the highest input min-entropy difference between OBC and NBC is code [{sorted_list_codes_by_OBC_NBC_difference[0][0]}, {sorted_list_codes_by_OBC_NBC_difference[0][1]}, {sorted_list_codes_by_OBC_NBC_difference[0][2]}] with difference of {highest_ent_diff_str}')
        f_hmin_diff.write(f'\nCode with the lowest input min-entropy difference between OBC and NBC is code [{sorted_list_codes_by_OBC_NBC_difference[-1][0]}, {sorted_list_codes_by_OBC_NBC_difference[-1][1]}, {sorted_list_codes_by_OBC_NBC_difference[-1][2]}] with difference of {lowest_ent_diff_str}')
        f_hmin_diff.write(f'\nCode with the highest input min-entropy relative reduction between OBC and NBC is code [{sorted_list_codes_by_OBC_NBC_relative_difference[0][0]}, {sorted_list_codes_by_OBC_NBC_relative_difference[0][1]}, {sorted_list_codes_by_OBC_NBC_relative_difference[0][2]}] with reduction of {sorted_list_codes_by_OBC_NBC_relative_difference[0][5]*100:.2f} %')
        f_hmin_diff.write(f'\nCode with the lowest input min-entropy relative reduction between OBC and NBC is code [{sorted_list_codes_by_OBC_NBC_relative_difference[-1][0]}, {sorted_list_codes_by_OBC_NBC_relative_difference[-1][1]}, {sorted_list_codes_by_OBC_NBC_relative_difference[-1][2]}] with reduction of {sorted_list_codes_by_OBC_NBC_relative_difference[-1][5]*100:.2f} %')
    return sorted_list_codes_by_OBC_NBC_difference, sorted_list_codes_by_OBC_NBC_relative_difference

# %%
'''
Applicable to: All codes
What: Determine which correctors are optimal
'''

def determine_Optimal_correctors(Entropy_bound_code_list, file_out_write_optimal, set_name):
    unravel_entropy_bound_code_list =  [list(t) for t in zip(*Entropy_bound_code_list)]

    Optimal_code_rate, Optimal_in_min_entropy, Optimal_n, Optimal_k, Optimal_d, Optimal_eff, Optimal_all = Optimal_frontier_2d_function(unravel_entropy_bound_code_list[3], unravel_entropy_bound_code_list[4], unravel_entropy_bound_code_list[0], unravel_entropy_bound_code_list[1], unravel_entropy_bound_code_list[2], unravel_entropy_bound_code_list[5], maxX = True, maxY = False)
    print(f'{set_name}: There are {len(Optimal_all)} optimal correctors.')

    with open(file_out_write_optimal, 'w') as file_iterator:
        for i in range(len(Optimal_all)):
            H_in_opt_str = str(Optimal_in_min_entropy[i].quantize(Decimal('.00000001'), rounding = ROUND_UP))
            efficiency_opt_str = str(Optimal_eff[i].quantize(Decimal('.00000001'), rounding = ROUND_HALF_EVEN))
            file_iterator.write(f'n = {Optimal_n[i]}, k = {Optimal_k[i]}, d = {Optimal_d[i]}, CR = {Optimal_code_rate[i]:.8f}, H_in_req = {H_in_opt_str}, efficiency = {efficiency_opt_str} \n')
        file_iterator.write(f'{set_name}: There are {len(Optimal_all)} optimal correctors.')
    return Optimal_code_rate, Optimal_in_min_entropy, Optimal_all

# %%
'''
Applicable to: All codes
Bounds: Lacharme and New bound 
What: Find Optimal codes for 9 given input min-entropies
'''
def determine_Optimal_correctors_for_targeted_h_in (Optimal_all, file_out_write_optimal_for_targeted_h_in, optimal_correctors_for_targeted_h_in_list, set_name, target_in_min_entropy_list,code_weight_distribution_list=[]):

    sorted_by_Hmin_Optimal_all = sorted(Optimal_all, key=lambda x: x[1], reverse=True) # sort by min-entropy

    with open(file_out_write_optimal_for_targeted_h_in, 'w') as fr:
        for target_in_min_entropy in target_in_min_entropy_list:
            found = 0
            for code in sorted_by_Hmin_Optimal_all:
                if code[1] <= target_in_min_entropy:
                    if 'NBC' in set_name:
                        for i in range(len(code_weight_distribution_list)):
                            if [code[2], code[3], code[4]] == code_weight_distribution_list[i][0]:
                                found_code_weight_distribution = code_weight_distribution_list[i][1]
                                break
                        curr_eff = NEW_H_min_bit_bound (code[2], code[3], Decimal(str(target_in_min_entropy)), found_code_weight_distribution) / (Decimal(str(target_in_min_entropy)) * code[2])
                    else:
                        curr_eff = OLD_H_min_bit_bound(code[3], Decimal(str(target_in_min_entropy)), code[4])/(Decimal(str(target_in_min_entropy)) * code[2])
                    curr_eff_str = str(curr_eff.quantize(Decimal('.00000001'), rounding = ROUND_HALF_EVEN)) # Regular rounding
                    H_in_req_str = str(code[1].quantize(Decimal('.00000001'), rounding = ROUND_UP))
                    fr.write(f'n = {code[2]}, k = {code[3]}, d = {code[4]}, CR = {code[0]:.8f}, H_in_req = {H_in_req_str}, H_in_target = {target_in_min_entropy}, targeted efficiency = {curr_eff_str} \n')
                    optimal_correctors_for_targeted_h_in_list.append(tuple(list(code) + [curr_eff]))
                    found = 1
                    break
            if found == 0:
                fr.write('WARNING: Adequate corrector not found!\n')

 ###################################################################################################################
 ###################################################################################################################
 ###################################################################################################################
 ###################################################################################################################
 ###################################################################################################################
 ###################################################################################################################
                

# sorted_by_Hmin_Pareto_NBCCYC_3D = sorted(Pareto_NBCCYC_3D, key=lambda x: x[4], reverse=True) # sort by min-entropy
# best_NBCCYC_3D_codes_list = list()
# best_NBCCYC_3D_codes_eff_list = list()

# file_write_targeted_Pareto_NBCCYC_3D = f'v2_Optimal_correctors_NBCCYC_3D_9_targeted_H_in_H_out_0.999.txt'

# with open(file_write_targeted_Pareto_NBCCYC_3D, 'w') as fr:
#     for target_in_min_entropy in target_in_min_entropy_list:
#         found = 0
#         for code in sorted_by_Hmin_Pareto_NBCCYC_3D:
#             if code[4] <= target_in_min_entropy:
#                 for i in range(len(code_weight_distribution_list)):
#                     if [code[0], code[1], code[2]] == code_weight_distribution_list[i][0]:
#                         found_code_weight_distribution = code_weight_distribution_list[i][1]
#                         break
#                 curr_eff = new_exact_H_min_bit_bound (code[0], code[1], target_in_min_entropy, found_code_weight_distribution)/(target_in_min_entropy * code[0])
#                 fr.write(f'n = {code[0]}, k = {code[1]}, d = {code[2]}, CR = {code[3]}, H_in_req = {np.ceil(code[4]*10**6)/(10**6)}, targeted efficiency = {curr_eff:.8f}, area = {code[7]:.2f} \n')
#                 best_NBCCYC_3D_codes_list.append(code)
#                 best_NBCCYC_3D_codes_eff_list.append(curr_eff)
#                 found = 1
#                 break
#         if found == 0:
#             fr.write('WARNING: Adequate code not found! Using simple XOR!\n')


# %%
def PLOT_abs_and_rel_entropy_differences_OBC_and_NBC (sorted_list_codes_by_OBC_NBC_difference, sorted_list_codes_by_OBC_NBC_relative_difference, code_weight_distribution_list, pdf_title):
    '''
    Find correctors with the biggest absolute and relative input min-entropy difference between the OBC and the NBC sets
    '''
    in_H_min_points = np.linspace(0, 1, num = 1000)

    n_abs_max_diff = sorted_list_codes_by_OBC_NBC_difference[0][0]
    k_abs_max_diff = sorted_list_codes_by_OBC_NBC_difference[0][1]
    d_abs_max_diff = sorted_list_codes_by_OBC_NBC_difference[0][2]

    for i in range(len(code_weight_distribution_list)):
        if [n_abs_max_diff, k_abs_max_diff, d_abs_max_diff] == code_weight_distribution_list[i][0]:
            code_weight_distribution_abs_max = code_weight_distribution_list[i][1]
            break

    OBC_abs_max_diff_vals_wo_bounding = [OLD_H_min_bit_bound (k_abs_max_diff, Decimal(str(in_H_min)), d_abs_max_diff) - k_abs_max_diff + 1 for in_H_min in in_H_min_points]
    OBC_abs_max_diff_vals = [max(0, float(str(OBC_abs_max_diff_vals_wo_bounding_i.quantize(Decimal('.00000000000000001'), rounding = ROUND_HALF_EVEN)))) for OBC_abs_max_diff_vals_wo_bounding_i in OBC_abs_max_diff_vals_wo_bounding]
    NBC_abs_max_diff_vals_wo_bounding = [NEW_H_min_bit_bound (n_abs_max_diff, k_abs_max_diff, Decimal(str(in_H_min)), code_weight_distribution_abs_max) - k_abs_max_diff + 1 for in_H_min in in_H_min_points]
    NBC_abs_max_diff_vals = [max(0, float(str(NBC_abs_max_diff_vals_wo_bounding_i.quantize(Decimal('.00000000000000001'), rounding = ROUND_HALF_EVEN)))) for NBC_abs_max_diff_vals_wo_bounding_i in NBC_abs_max_diff_vals_wo_bounding]


    n_rel_max_diff = sorted_list_codes_by_OBC_NBC_relative_difference[0][0]
    k_rel_max_diff = sorted_list_codes_by_OBC_NBC_relative_difference[0][1]
    d_rel_max_diff = sorted_list_codes_by_OBC_NBC_relative_difference[0][2]


    for i in range(len(code_weight_distribution_list)):
        if [n_rel_max_diff, k_rel_max_diff, d_rel_max_diff] == code_weight_distribution_list[i][0]:
            code_weight_distribution_rel_max = code_weight_distribution_list[i][1]
            break


    OBC_rel_max_diff_vals_wo_bounding = [OLD_H_min_bit_bound (k_rel_max_diff, Decimal(str(in_H_min)), d_rel_max_diff) - k_rel_max_diff + 1 for in_H_min in in_H_min_points]
    OBC_rel_max_diff_vals = [max(0, float(str(OBC_rel_max_diff_vals_wo_bounding_i.quantize(Decimal('.00000000000000001'), rounding = ROUND_HALF_EVEN)))) for OBC_rel_max_diff_vals_wo_bounding_i in OBC_rel_max_diff_vals_wo_bounding]
    NBC_rel_max_diff_vals_wo_bounding = [NEW_H_min_bit_bound (n_rel_max_diff, k_rel_max_diff, Decimal(str(in_H_min)), code_weight_distribution_rel_max) - k_rel_max_diff + 1 for in_H_min in in_H_min_points]
    NBC_rel_max_diff_vals = [max(0, float(str(NBC_rel_max_diff_vals_wo_bounding_i.quantize(Decimal('.00000000000000001'), rounding = ROUND_HALF_EVEN)))) for NBC_rel_max_diff_vals_wo_bounding_i in NBC_rel_max_diff_vals_wo_bounding]


    '''
    Draw figure with output-input min-entropy for both bounds and
    codes with the biggest min-entropy differences
    '''
    with plt.style.context(['science']):

        width_pt = ieee_one_column_inch*72 + 36 #264#307.34963919#264 # Adjust this number so that the output image has width of 3.5 in

        legend_fontsize = 8
        title_size = 10
        axes_label_size = 8
        axes_tick_size = 8
        scatter_marker_size = 2
        makerscale_factor = 3

        plt.rc('xtick', labelsize = axes_tick_size)
        plt.rc('ytick', labelsize = axes_tick_size)
        plt.rc('axes', labelsize = axes_label_size)
        plt.rc('legend', fontsize = legend_fontsize)
        

        fig = plt.figure(figsize = set_size(width_pt), dpi = 600)
        text_fontsize = 7

        l1, = plt.plot(in_H_min_points, NBC_rel_max_diff_vals, color = '#cc0000', label = f'Code [{n_rel_max_diff}, {k_rel_max_diff}, {d_rel_max_diff}]:\nNew bound', zorder=1)
        l2, = plt.plot(in_H_min_points, OBC_rel_max_diff_vals, color = '#00cccc' , label = f'Code [{n_rel_max_diff}, {k_rel_max_diff}, {d_rel_max_diff}]:\nOld bound', zorder=2, linestyle='--')
        m1 = plt.scatter([0.27445],[0.999],c = 'black', marker = 'x', s = 10, zorder=3)
        plt.text(0.09, 0.93, r'$\left(0.275, 0.999\right)$', fontsize = text_fontsize, color = 'black', bbox=dict(facecolor='white', edgecolor='none', pad=0))
        m2 = plt.scatter([0.71507],[0.999],c = 'black', marker = '^', s = 10, zorder=5)
        plt.text(0.56, 0.93, r'$\left(0.715, 0.999\right)$', fontsize = text_fontsize, color = 'black', bbox=dict(facecolor='white', edgecolor='none', pad=0))
        

        l3, = plt.plot(in_H_min_points, NBC_abs_max_diff_vals, color = '#66cc00', label = f'Code [{n_abs_max_diff}, {k_abs_max_diff}, {d_abs_max_diff}]:\nNew bound', zorder=1, linestyle=':')
        l4, = plt.plot(in_H_min_points, OBC_abs_max_diff_vals, color = '#6600cc' , label = f'Code [{n_abs_max_diff}, {k_abs_max_diff}, {d_abs_max_diff}]:\nOld bound', zorder=2, linestyle='-.')
        m3 = plt.scatter([0.408],[0.999],c = 'black', marker = '*', s = 10, zorder=3)
        plt.text(0.32, 0.93, r'$\left(0.408, 0.999\right)$', fontsize = text_fontsize, color = 'black', bbox=dict(facecolor='white', edgecolor='none', pad=0))
        m4 = plt.scatter([0.8543],[0.999],c = 'black', marker = 'o', s = 10, zorder=5)
        plt.text(0.79, 0.93, r'$\left(0.854, 0.999\right)$', fontsize = text_fontsize, color = 'black', bbox=dict(facecolor='white', edgecolor='none', pad=0))


        plt.xlabel(r'$\mathrm{H}^{in}_{\infty}$')
        plt.ylabel(r'$\mathrm{H}^{out,\, 1}_{\infty}$')
    
        plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))

        handles=[(l1,m1),(l2,m2),(l3,m3),(l4,m4)]
        _, labels = fig.get_axes()[0].get_legend_handles_labels()

        plt.legend(bbox_to_anchor =(0.5,-0.55), loc='lower center', ncol=2, fontsize = legend_fontsize, labels = labels, handles=handles, handler_map = {tuple: matplotlib.legend_handler.HandlerTuple(None)}, handlelength=3)

        plt.savefig(pdf_title, format='pdf', dpi = 600, bbox_inches='tight', pad_inches=0.005)


# %%
def find_efficiencies_for_entire_h_in_range (Optimal_OBC_correctors, Optimal_NBC_correctors, code_weight_distribution_list):

    sorted_by_Hmin_Optimal_OBC = sorted(Optimal_OBC_correctors, key=lambda x: x[1], reverse=True) # sort by min-entropy
    added_min_entropy_list = list()
    added_efficiency_list = list()

    # Make list of equally spaced points between 0 and the minimum output entropy per bit bound
    target_in_min_entropy_list = np.linspace(0, float(minimum_output_entropy_per_bit_bound), 100001)

    for target_in_min_entropy in target_in_min_entropy_list:
        for code in sorted_by_Hmin_Optimal_OBC:
            if code[1] <= target_in_min_entropy:
                curr_eff = OLD_H_min_bit_bound(code[3], Decimal(str(target_in_min_entropy)), code[4])/(Decimal(str(target_in_min_entropy)) * code[2])
                added_min_entropy_list.append(target_in_min_entropy)
                added_efficiency_list.append(curr_eff)
                break


    OBC_in_min_entropies_range = added_min_entropy_list
    Optimal_OBC_efficiencies = added_efficiency_list


    sorted_by_Hmin_Pareto_NBC = sorted(Optimal_NBC_correctors, key=lambda x: x[1], reverse=True) # sort by min-entropy
    added_min_entropy_list = list()
    added_efficiency_list = list()


    # for target_in_min_entropy in target_in_min_entropy_list:
    for target_in_min_entropy in OBC_in_min_entropies_range:
        for code in sorted_by_Hmin_Pareto_NBC:
            if code[1] <= target_in_min_entropy:
                for i in range(len(code_weight_distribution_list)):
                    if [code[2], code[3], code[4]] == code_weight_distribution_list[i][0]:
                        found_code_weight_distribution = code_weight_distribution_list[i][1]
                        break
                curr_eff = NEW_H_min_bit_bound (code[2], code[3], Decimal(str(target_in_min_entropy)), found_code_weight_distribution)/(Decimal(str(target_in_min_entropy)) * code[2])                    
                added_min_entropy_list.append(target_in_min_entropy)
                added_efficiency_list.append(curr_eff)
                break

    NBC_in_min_entropies_range =  added_min_entropy_list
    Optimal_NBC_efficiencies = added_efficiency_list

    return OBC_in_min_entropies_range, Optimal_OBC_efficiencies, NBC_in_min_entropies_range, Optimal_NBC_efficiencies


# %%
def PLOT_OBC_NBC_correctors_and_efficiencies (OBC_entropy_bound_code_list, NBC_entropy_bound_code_list, Optimal_OBC_correctors_in_min_entropy, Optimal_NBC_correctors_in_min_entropy, Optimal_OBC_correctors_code_rate,  Optimal_NBC_correctors_code_rate, OBC_in_min_entropies_range, Optimal_OBC_efficiencies, NBC_in_min_entropies_range, Optimal_NBC_efficiencies, pdf_title):
    '''Generate plots for the paper'''
    with plt.style.context(['science']):
        width_pt = ieee_two_column_inch*73 + 16 #307.34963919#524
        legend_fontsize = 8
        title_size = 10
        axes_label_size = 8
        scatter_marker_size = 2
        makerscale_factor = 3

        plt.rc('xtick', labelsize = 8)
        plt.rc('ytick', labelsize = 8)
        plt.rc('axes', labelsize = axes_label_size)
        plt.rc('legend', fontsize = legend_fontsize)
        

        dim_0, dim_1 = set_size(width_pt, subplots = (2,2))

        fig = plt.figure(figsize=(np.round(dim_0,2),np.round(dim_1,2)), dpi = 600)
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)

        ax1.set_title('(a)', size = axes_label_size)
        ax1.scatter([i[-2] for i in OBC_entropy_bound_code_list], [i[-3] for i in OBC_entropy_bound_code_list], c = '#00cccc', marker = 'o', s = scatter_marker_size, label = 'Appropriate correctors from OBC')
        ax1.plot(Optimal_OBC_correctors_in_min_entropy, [i/minimum_output_entropy_per_bit_bound for i in Optimal_OBC_correctors_in_min_entropy], 'k-.', label = 'Extraction limit')
        ax1.set_ylim((-0.02,1.02))
        ax1.set_xlim((-0.02,1.02))
        ax1.set_xlabel(r'$\mathrm{H}^{in,\, req}_{\infty}$')
        ax1.set_ylabel(r'Code rate')
        ax1.legend(loc = 'upper left', markerscale = makerscale_factor, fontsize = legend_fontsize)

        ax2.set_title('(b)', size = axes_label_size)
        ax2.scatter([i[-2] for i in NBC_entropy_bound_code_list], [i[-3] for i in NBC_entropy_bound_code_list], c = '#cc0000', marker = 'o', s = scatter_marker_size, label = 'Appropriate correctors from NBC')
        ax2.plot(Optimal_NBC_correctors_in_min_entropy, [i/minimum_output_entropy_per_bit_bound for i in Optimal_NBC_correctors_in_min_entropy], 'k-.', label = 'Extraction limit')
        ax2.set_ylim((-0.02,1.02))
        ax2.set_xlim((-0.02,1.02))
        ax2.set_xlabel(r'$\mathrm{H}^{in,\, req}_{\infty}$')
        ax2.set_ylabel(r'Code rate')
        ax2.legend(loc = 'upper left', markerscale = makerscale_factor, fontsize = legend_fontsize)

        ax3.set_title('(c)', size = axes_label_size)
        ax3.scatter(Optimal_OBC_correctors_in_min_entropy, Optimal_OBC_correctors_code_rate, c = '#00cccc', marker = 's', s = scatter_marker_size, label = f'OBC PF correctors ({len(Optimal_OBC_correctors_in_min_entropy)})')
        ax3.scatter(Optimal_NBC_correctors_in_min_entropy, Optimal_NBC_correctors_code_rate, c = '#cc0000', marker = 's', s = scatter_marker_size, label = f'NBC PF correctors ({len(Optimal_NBC_correctors_in_min_entropy)})')
        ax3.plot(Optimal_NBC_correctors_in_min_entropy, [i/minimum_output_entropy_per_bit_bound for i in Optimal_NBC_correctors_in_min_entropy], 'k-.', label = 'Extraction limit')
        ax3.set_ylim((-0.02,1.02))
        ax3.set_xlim((-0.02,1.02))
        ax3.set_xlabel(r'$\mathrm{H}^{in,\, req}_{\infty}$')
        ax3.set_ylabel(r'Code rate')
        ax3.legend(loc = 'upper left', markerscale = makerscale_factor, fontsize = legend_fontsize)

        ax4.set_title('(d)', size = axes_label_size)

        l1 = ax4.plot(OBC_in_min_entropies_range, Optimal_OBC_efficiencies, c = '#00cccc', label = f'Old bound', zorder=1)
        l2 = ax4.plot(NBC_in_min_entropies_range, Optimal_NBC_efficiencies, c = '#cc0000', label = f'New bound', zorder=2)
        
        axins = ax4.inset_axes([0.15,0.05,0.42,0.45])
        axins.plot(OBC_in_min_entropies_range, Optimal_OBC_efficiencies, c = '#00cccc')
        axins.plot(NBC_in_min_entropies_range, Optimal_NBC_efficiencies, c = '#cc0000')
        # subregion of the original image
        x1, x2, y1, y2 = 0.03946, 0.039495, 0.18, 0.57
        mark_inset(ax4, axins, loc1=1, loc2=3, fc="none", ec= 'black', lw = 0.5)
        
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_xticklabels([])
        axins.set_yticklabels([])


        ax4.set_ylim((-0.02,1.02))
        ax4.set_xlim((-0.02,1.02))

        ax4.set_xlabel(r'$\mathrm{H}^{in}_{\infty}$')
        ax4.set_ylabel(r'Extraction efficiency $\left(\eta\right)$')


        m1 = ax4.scatter([0.76697226],[0.5760907909871946],c = 'black', marker = "^", s = 10, zorder=3)
        ax4.text(0.6, 0.48, r'$\left(0.76697, 0.5761\right)$', fontsize = 7, color = 'black')
        ax4.scatter([0.76697226],[0.9727754695905784],c = 'black', marker = "x", s = 10, zorder=4)
        m2 = ax4.text(0.47, 0.94, r'$\left(0.76697, 0.9728\right)$', fontsize = 7, color = 'black')

        axins.scatter([0.03947049],[0.2146046523342227],c = 'black', marker = "^", s = 10, zorder=3)
        axins.annotate(r'$\left(0.03947, 0.2146\right)$', xy=(0.039471, 0.21), fontsize=7, color = 'black')
        axins.scatter([0.03947049],[0.49478229043530364],c = 'black', marker = "x", s = 10, zorder=4)
        axins.annotate(r'$\left(0.03947, 0.4948\right)$', xy=(0.039471, 0.52), fontsize=7, color = 'black')

        handles=[(l1,m1),(l2,m2)]
        _, labels = ax4.get_legend_handles_labels()

        ax4.legend(loc = 'lower right', markerscale = makerscale_factor, fontsize = legend_fontsize)

        plt.tight_layout()
        plt.savefig(pdf_title, format='pdf', dpi = 600, bbox_inches='tight', pad_inches=0.005)   
    


# %%
def PLOT_OBCCYC_NBCCYC_correctors (Optimal_OBCCYC_correctors_in_min_entropy, Optimal_NBCCYC_correctors_in_min_entropy, Optimal_OBCCYC_correctors_code_rate,  Optimal_NBCCYC_correctors_code_rate, pdf_title):
    '''
    Draw figure with output-input min-entropy for both bounds and
    codes with the biggest min-entropy differences
    '''
    with plt.style.context(['science']):
        width_pt = ieee_one_column_inch*72 + 39

        legend_fontsize = 8
        title_size = 10
        axes_label_size = 8
        scatter_marker_size = 2
        makerscale_factor = 3

        plt.rc('xtick', labelsize = 8)
        plt.rc('ytick', labelsize = 8)
        plt.rc('axes', labelsize = axes_label_size)
        plt.rc('legend', fontsize = legend_fontsize)


        fig = plt.figure(figsize = set_size(width_pt), dpi = 600)
        text_fontsize = 6


        plt.scatter(Optimal_OBCCYC_correctors_in_min_entropy, Optimal_OBCCYC_correctors_code_rate, c = '#00cccc', marker = 's', s = scatter_marker_size, label = f'OBCCYC PF correctors ({len(Optimal_OBCCYC_correctors_in_min_entropy)})')
        plt.scatter(Optimal_NBCCYC_correctors_in_min_entropy, Optimal_NBCCYC_correctors_code_rate, c = '#cc0000', marker = 's', s = scatter_marker_size, label = f'NBCCYC PF correctors ({len(Optimal_NBCCYC_correctors_in_min_entropy)})')
        plt.plot(Optimal_NBCCYC_correctors_in_min_entropy, [i/minimum_output_entropy_per_bit_bound for i in Optimal_NBCCYC_correctors_in_min_entropy], 'k-.', label = 'Extraction limit')
        plt.ylim((-0.02,1.02))
        plt.xlim((-0.02,1.02))
        plt.xlabel(r'$\mathrm{H}^{in,\, req}_{\infty}$')
        plt.ylabel(r'Code rate')
        plt.legend(loc = 'upper left', markerscale = makerscale_factor, fontsize = legend_fontsize)

        plt.savefig(pdf_title, format='pdf', dpi = 600, bbox_inches='tight', pad_inches=0.005)


# %%
def PLOT_OBCCYC_NBCCYC_efficiencies (OBCCYC_in_min_entropies_range, Optimal_OBCCYC_efficiencies, NBCCYC_in_min_entropies_range, Optimal_NBCCYC_efficiencies, pdf_title):
    '''
    Draw figure with output-input min-entropy for both bounds and
    codes with the biggest min-entropy differences
    '''
    with plt.style.context(['science']):
        width_pt = ieee_one_column_inch*72 + 37.8

        legend_fontsize = 8
        title_size = 10
        axes_label_size = 8
        scatter_marker_size = 2
        makerscale_factor = 3
        text_on_graph_size = 7

        plt.rc('xtick', labelsize = 8)
        plt.rc('ytick', labelsize = 8)
        plt.rc('axes', labelsize = axes_label_size)
        plt.rc('legend', fontsize = legend_fontsize)


        dim_0, dim_1 = set_size(width_pt, subplots = (2,2))
        fig = plt.figure(figsize=(np.round(dim_0,2),np.round(dim_1+0.017,3)), dpi = 600)

        #fig = plt.figure(figsize = set_size(width_pt), dpi = 600)
        ax4 = fig.add_subplot(111)
        text_fontsize = 6


        ax4.plot(OBCCYC_in_min_entropies_range, Optimal_OBCCYC_efficiencies, c = '#00cccc', label = f'Old bound', zorder=1)
        ax4.plot(NBCCYC_in_min_entropies_range, Optimal_NBCCYC_efficiencies, c = '#cc0000', label = f'New bound', zorder=2)
        
        axins = ax4.inset_axes([0.15,0.05,0.42,0.45])
        axins.plot(OBCCYC_in_min_entropies_range, Optimal_OBCCYC_efficiencies, c = '#00cccc')
        axins.plot(NBCCYC_in_min_entropies_range, Optimal_NBCCYC_efficiencies, c = '#cc0000')
        # subregion of the original image
        x1, x2, y1, y2 = 0.039545, 0.039565, 0.13, 0.53
        mark_inset(ax4, axins, loc1=1, loc2=3, fc="none", ec= 'black', lw = 0.5)
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_xticklabels([])
        axins.set_yticklabels([])
        axins.scatter([0.039550409999999994],[0.17541585483096894],c = 'black', marker = "^", s = 10, zorder=3)
        axins.annotate(r'$\left(0.03955, 0.1754\right)$', xy=(0.03955, 0.21), fontsize=text_on_graph_size, color = 'black')
        axins.scatter([0.039550409999999994],[0.49474883266909225],c = 'black', marker = "x", s = 10, zorder=4)
        axins.annotate(r'$\left(0.03955, 0.4947\right)$', xy=(0.03955, 0.45), fontsize=text_on_graph_size, color = 'black')


        ax4.set_ylim((-0.02,1.02))
        ax4.set_xlim((-0.02,1.02))
        ax4.set_xlabel(r'$\mathrm{H}^{in}_{\infty}$')
        ax4.set_ylabel(r'Extraction efficiency $\left(\eta\right)$')
        ax4.legend(loc = 'lower right', markerscale = makerscale_factor, fontsize = legend_fontsize)

        ax4.scatter([0.8843048099999999],[0.5899909798881077],c = 'black', marker = "^", s = 10, zorder=3)
        ax4.text(0.7, 0.48, r'$\left(0.8843, 0.59\right)$', fontsize = text_on_graph_size, color = 'black')
        ax4.scatter([0.8843048099999999],[0.9534419880646975],c = 'black', marker = "x", s = 10, zorder=4)
        ax4.text(0.56, 0.94, r'$\left(0.8843, 0.9534\right)$', fontsize = text_on_graph_size, color = 'black')
        
        plt.savefig(pdf_title, format='pdf', dpi = 600, bbox_inches='tight', pad_inches=0.005)