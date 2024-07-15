def check_early_stop_condition(inputs, logfile, epoch, train_epoch_result, valid_epoch_result, early_stop_dict):
    if 'early_stop_criteria' not in inputs['neural_network']:
        return False

    # Save RMSE values for all convergence conditions
    update_early_stop_dict(inputs, epoch, train_epoch_result, valid_epoch_result, early_stop_dict)

    criteria_inputs = inputs['neural_network']['early_stop_criteria']
    avg_range = criteria_inputs['averaging']
    patience = criteria_inputs['patience']

    # Pass until minimum epoch training done
    if epoch < patience + avg_range:
        return False

    # Check all convergence conditions are satisfied
    if len(early_stop_dict.keys()) != 0:
        pass_key = 1
    else:
        pass_key = 0

    for key in early_stop_dict.keys():
        l = len(early_stop_dict[key][0])
        current_rmse = sum(early_stop_dict[key][0][l-avg_range:])/avg_range
        prev_rmse = sum(early_stop_dict[key][0][l-avg_range-patience:l-patience])/avg_range

        if abs(prev_rmse - current_rmse) > early_stop_dict[key][1]:
            pass_key = 0
            break

    if pass_key == 1:
        return True
    else:
        return False

def update_early_stop_dict(inputs, epoch, train_epoch_result, valid_epoch_result, early_stop_dict):
    criteria_inputs = inputs['neural_network']['early_stop_criteria']

    if 'train_set_condition' in criteria_inputs:
        epoch_result = train_epoch_result
        str_labels = criteria_inputs['train_set_condition']['target_structures']
        E_threshold = criteria_inputs['train_set_condition']['E_threshold'] if 'E_threshold' in criteria_inputs['train_set_condition'] else None
        F_threshold = criteria_inputs['train_set_condition']['F_threshold'] if 'F_threshold' in criteria_inputs['train_set_condition'] else None
        S_threshold = criteria_inputs['train_set_condition']['S_threshold'] if 'S_threshold' in criteria_inputs['train_set_condition'] else None
        for str_label in str_labels:
            if str_label == 'Total':
                if E_threshold:
                    append_avg_rmse(early_stop_dict, epoch_result, E_threshold, key='e_err', label='train_e')
                if F_threshold:
                    append_avg_rmse(early_stop_dict, epoch_result, F_threshold, key='f_err', label='train_f')
                if S_threshold:
                    append_avg_rmse(early_stop_dict, epoch_result, S_threshold, key='s_err', label='train_s')
            else:
                if str_label not in epoch_result['struct_labels']:
                    continue

                if E_threshold:
                    append_structure_rmse(early_stop_dict, epoch_result, E_threshold, key='e_err', str_label=str_label, label='train_e_'+str_label)
                if F_threshold:
                    append_structure_rmse(early_stop_dict, epoch_result, F_threshold, key='f_err', str_label=str_label, label='train_f_'+str_label)
                if S_threshold:
                    append_structure_rmse(early_stop_dict, epoch_result, S_threshold, key='s_err', str_label=str_label, label='train_s_'+str_label)

    if 'valid_set_condition' in criteria_inputs:
        epoch_result = valid_epoch_result
        str_labels = criteria_inputs['valid_set_condition']['target_structures']
        E_threshold = criteria_inputs['valid_set_condition']['E_threshold'] if 'E_threshold' in criteria_inputs['valid_set_condition'] else None
        F_threshold = criteria_inputs['valid_set_condition']['F_threshold'] if 'F_threshold' in criteria_inputs['valid_set_condition'] else None
        S_threshold = criteria_inputs['valid_set_condition']['S_threshold'] if 'S_threshold' in criteria_inputs['valid_set_condition'] else None
        for str_label in str_labels:
            if str_label == 'Total':
                if E_threshold:
                    append_avg_rmse(early_stop_dict, epoch_result, E_threshold,key='e_err', label='valid_e')
                if F_threshold:
                    append_avg_rmse(early_stop_dict, epoch_result, F_threshold, key='f_err', label='valid_f')
                if S_threshold:
                    append_avg_rmse(early_stop_dict, epoch_result, S_threshold, key='s_err', label='valid_s')
            else:
                if str_label not in epoch_result['struct_labels']:
                    continue

                if E_threshold:
                    append_structure_rmse(early_stop_dict, epoch_result, E_threshold, key='e_err', str_label=str_label, label='valid_e_'+str_label)
                if F_threshold:
                    append_structure_rmse(early_stop_dict, epoch_result, F_threshold, key='f_err', str_label=str_label, label='valid_f_'+str_label)
                if S_threshold:
                    append_structure_rmse(early_stop_dict, epoch_result, S_threshold, key='s_err', str_label=str_label, label='valid_s_'+str_label)

def append_avg_rmse(early_stop_dict, epoch_result, stop_value, key='e_err', label=None):
    if label not in early_stop_dict.keys():
        early_stop_dict[label] = ([], stop_value)

    v_sum = 0
    v_count = 0
    for str_label in epoch_result['struct_labels']:
        v_sum   += epoch_result[key][str_label].sum
        v_count += epoch_result[key][str_label].count
    v_rmse = (v_sum / v_count) ** 0.5

    early_stop_dict[label][0].append(v_rmse)

def append_structure_rmse(early_stop_dict, epoch_result, stop_value, key='e_err', str_label=None, label=None):
    if label not in early_stop_dict.keys():
        early_stop_dict[label] = ([], stop_value)

    v_sum = epoch_result[key][str_label].sum
    v_count = epoch_result[key][str_label].count
    v_rmse = (v_sum / v_count) ** 0.5

    early_stop_dict[label][0].append(v_rmse)

