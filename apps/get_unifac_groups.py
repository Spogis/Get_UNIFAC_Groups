import operator
from apps.fragmenter import fragmenter
from apps.extract_data_from_log import *
from apps.smarts import UNIFAC_SMARTS
from rdkit import Chem



def Get_UNIFAC_Groups():
    my_dataset = 'apps/reference_DB.csv'

    # def info_to_CSV(inchikey, SMILES, pubchem_id, fragmentation):
    #     fragmentation_array = []
    #     for group_number, amount in fragmentation.items():
    #         fragmentation_array.append(str(group_number) + ":" + str(amount))
    #
    #     return inchikey + "," + SMILES + "," + pubchem_id + "," + "|".join(fragmentation_array)
    #
    def CSV_to_info(CSV_line, has_fragmentation=False):
        CSV_line = CSV_line.replace('\n', '')
        array = CSV_line.split(',')

        fragmentation = {}

        if has_fragmentation:
            fragmentation_array = array[3].split('|')
            for match_str in fragmentation_array:
                array2 = match_str.split(':')
                group_number = int(array2[0])
                amount = int(array2[1])

                fragmentation[group_number] = amount

        return array[0], array[1], array[2], fragmentation

    def function_to_choose_fragmentation(fragmentations):
        fragmentations_descriptors = {}
        i = 0
        for fragmentation in fragmentations:
            fragmentations_descriptors[i] = [len(fragmentation)]
            i += 1

        sorted_fragmentations_dict = sorted(fragmentations_descriptors.items(), key=operator.itemgetter(1))

        return fragmentations[sorted_fragmentations_dict[0][0]]


    def is_fragmentation_equal_to_other_fragmentation(fragmentation, other_fragmentation):
        for group_number, amount in fragmentation.items():
            if group_number in other_fragmentation:
                if fragmentation[group_number] != other_fragmentation[group_number]:
                    return False
        return True

    def log_structure_results(f, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB,
                              status=''):
        f.write('https://pubchem.ncbi.nlm.nih.gov/compound/' + pubchem_id + '#section=2D-Structure\n')
        f.write(SMILES + '\n')
        f.write(inchikey + '\n')
        f.write('\n' + 'Fragmentation was successfull: ' + str(success) + '\n')

        if status != '':
            f.write(status + '\n')

        if success:
            f.write('Fragmentation from the algorithm:\n')
            sorted_group_number = sorted(fragmentation.keys())

            for group_number in sorted_group_number:
                f.write((UNIFAC_SMARTS[group_number - 1][0]).ljust(12, ' ') + '\t' + str(group_number).ljust(8, ' ') + str(
                    fragmentation[group_number]).ljust(8, ' ') + '\n')

        f.write('\n')

        if len(fragmentation_reference_DB) > 0:
            f.write('Fragmentation from the reference database:\n')
            sorted_group_number = sorted(fragmentation_reference_DB.keys())

            for group_number in sorted_group_number:
                f.write((UNIFAC_SMARTS[group_number - 1][0]).ljust(12, ' ') + '\t' + str(group_number).ljust(8, ' ') + str(
                    fragmentation_reference_DB[group_number]).ljust(8, ' ') + '\n')

        f.write('\n\n')

    # get the fragmentation scheme in the format necessary
    fragmentation_scheme = {i + 1: j[1] for i, j in enumerate(UNIFAC_SMARTS)}

    # sort the fragmentation scheme according to the descriptors
    pattern_descriptors = {}
    for group_number, SMARTS in fragmentation_scheme.items():
        if type(SMARTS) is list:
            SMARTS = SMARTS[0]

        if SMARTS != "":
            pattern = fragmenter.get_mol_with_properties_from_SMARTS(SMARTS)

            pattern_descriptors[group_number] = [pattern.GetUnsignedProp('n_available_bonds') == 0, \
                                                 (pattern.GetBoolProp('is_simple_atom_on_c') or pattern.GetBoolProp(
                                                     'is_simple_atom')), \
                                                 pattern.GetUnsignedProp('n_atoms_defining_SMARTS'),
                                                 pattern.GetUnsignedProp('n_available_bonds') == 1, \
                                                 fragmenter.get_heavy_atom_count(pattern) - pattern.GetUnsignedProp(
                                                     'n_carbons'), \
                                                 pattern.GetBoolProp('has_atoms_in_ring'), \
                                                 pattern.GetUnsignedProp('n_triple_bonds'), \
                                                 pattern.GetUnsignedProp('n_double_bonds')]

    sorted_pattern_descriptors = sorted(pattern_descriptors.items(), key=operator.itemgetter(1), reverse=True)
    sorted_group_numbers = [i[0] for i in sorted_pattern_descriptors]

    # first step: fragment reference database and compare with the results
    reference_DB = []
    with open(my_dataset) as f:
        for line in f.readlines():
            reference_DB.append(CSV_to_info(line, True))

    reference_DB_fragmentation_stats = {}

    simple_fragmenter_fragmented = []
    simple_fragmenter_fragmented_and_equal_to_reference_DB = []
    complete_fragmenter_fragmented = []
    complete_fragmenter_fragmented_and_equal_to_reference_DB = []

    right_size_for_complete_fragmenter = []

    simple_fragmenter = fragmenter(fragmentation_scheme, 'simple')
    complete_fragmenter = fragmenter(fragmentation_scheme, 'complete', 100, function_to_choose_fragmentation)

    # without sorting the patterns
    # print('####################################################################')
    # print('Fragmenting the reference database without the patterns sorted (simple and complete algorithm)')

    i_structure = 0
    f_simple = open('logs/simple_fragmentation_without_patterns_sorted_results.log', 'w+')
    f_complete = open('logs/complete_fragmentation_without_patterns_sorted_results.log', 'w+')

    for inchikey, SMILES, pubchem_id, fragmentation_reference_DB in reference_DB:

        i_structure = i_structure + 1
        if i_structure % 2000 == 0:
            print('{:2.1f} .'.format((100.0 * i_structure) / len(reference_DB)), end=" ")

        lines = []

        for group_number, amount in fragmentation_reference_DB.items():
            if not group_number in reference_DB_fragmentation_stats:
                reference_DB_fragmentation_stats[group_number] = 0

            reference_DB_fragmentation_stats[group_number] += amount

        fragmentation, success = simple_fragmenter.fragment(SMILES)
        if success:
            simple_fragmenter_fragmented.append(inchikey)
            if is_fragmentation_equal_to_other_fragmentation(fragmentation, fragmentation_reference_DB):
                simple_fragmenter_fragmented_and_equal_to_reference_DB.append(inchikey)

        log_structure_results(f_simple, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB)

        n_heavy_atoms = 0
        for sub_SMILES in SMILES.split("."):
            n_heavy_atoms = max(n_heavy_atoms, fragmenter.get_heavy_atom_count(Chem.MolFromSmiles(sub_SMILES)))

        if n_heavy_atoms <= 20:
            right_size_for_complete_fragmenter.append(inchikey)
            fragmentation, success = complete_fragmenter.fragment(SMILES)
            if success:
                complete_fragmenter_fragmented.append(inchikey)
                if is_fragmentation_equal_to_other_fragmentation(fragmentation, fragmentation_reference_DB):
                    complete_fragmenter_fragmented_and_equal_to_reference_DB.append(inchikey)

            log_structure_results(f_complete, pubchem_id, SMILES, inchikey, success, fragmentation,
                                  fragmentation_reference_DB)
        else:
            log_structure_results(f_complete, pubchem_id, SMILES, inchikey, success, fragmentation,
                                  fragmentation_reference_DB, 'Structure was skipped because it is larger than 20 atoms.')

    f_simple.close()
    f_complete.close()

    right_size_for_complete_fragmenter2 = []
    # with sorting the patterns
    simple_fragmenter.fragmentation_scheme_order = sorted_group_numbers
    complete_fragmenter.fragmentation_scheme_order = sorted_group_numbers

    simple_fragmenter_sorted_fragmented = []
    simple_fragmenter_sorted_fragmented_and_equal_to_reference_DB = []
    complete_fragmenter_sorted_fragmented = []
    complete_fragmenter_sorted_fragmented_and_equal_to_reference_DB = []

    # print('####################################################################')
    # print('Fragmenting the reference database with the patterns sorted (simple and complete algorithm)')
    i_structure = 0
    f_simple = open('logs/simple_fragmentation_with_patterns_sorted_results.log', 'w+')
    f_complete = open('logs/complete_fragmentation_with_patterns_sorted_results.log', 'w+')
    for inchikey, SMILES, pubchem_id, fragmentation_reference_DB in reference_DB:

        i_structure = i_structure + 1

        if i_structure % 2000 == 0:
            print('{:2.1f} .'.format((100.0 * i_structure) / len(reference_DB)), end=" ")

        fragmentation, success = simple_fragmenter.fragment(SMILES)
        if success:
            simple_fragmenter_sorted_fragmented.append(inchikey)
            if is_fragmentation_equal_to_other_fragmentation(fragmentation, fragmentation_reference_DB):
                simple_fragmenter_sorted_fragmented_and_equal_to_reference_DB.append(inchikey)

        log_structure_results(f_simple, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB)

        n_heavy_atoms = 0
        for sub_SMILES in SMILES.split("."):
            n_heavy_atoms = max(n_heavy_atoms, fragmenter.get_heavy_atom_count(Chem.MolFromSmiles(sub_SMILES)))

        if n_heavy_atoms <= 20:
            right_size_for_complete_fragmenter2.append(inchikey)
            fragmentation, success = complete_fragmenter.fragment(SMILES)
            if success:
                complete_fragmenter_sorted_fragmented.append(inchikey)
                if is_fragmentation_equal_to_other_fragmentation(fragmentation, fragmentation_reference_DB):
                    complete_fragmenter_sorted_fragmented_and_equal_to_reference_DB.append(inchikey)

            log_structure_results(f_complete, pubchem_id, SMILES, inchikey, success, fragmentation,
                                  fragmentation_reference_DB)
        else:
            log_structure_results(f_complete, pubchem_id, SMILES, inchikey, success, fragmentation,
                                  fragmentation_reference_DB, 'Structure was skipped because it is larger than 20 atoms.')

    f_simple.close()
    f_complete.close()

    #log_file_path = 'logs/complete_fragmentation_with_patterns_sorted_results.log'
    log_file_path = 'logs/complete_fragmentation_without_patterns_sorted_results.log'
    #log_file_path = 'logs/simple_fragmentation_with_patterns_sorted_results.log'
    #log_file_path = 'logs/simple_fragmentation_without_patterns_sorted_results.log'

    return create_dataset_from_log(log_file_path)
