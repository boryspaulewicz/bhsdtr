// 'softmax' gamma link function
criteria = criteria_scale * inv_Phi(head(cumulative_sum(softmax(append_row(gamma, 0))), K - 1));
