//gamma_link:softmax
criteria = criteria_scale * inv_Phi(head(cumulative_sum(softmax(append_row(gamma, 0))), K - 1));
