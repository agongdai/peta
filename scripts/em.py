def fac(n):
    if n == 0:
        return 1
    else:
        return n * (fac(n-1))
    
def select(n, c):
    #print '%d, %d: %d' % (n, c, fac(n) / fac(c) / fac(n - c))
    return fac(n) / fac(c) / fac(n - c)

def dices():
    n_rounds = 5
    n_trials = 10
    n_heads = [5, 9, 8, 4, 7]
    theta_a = 0.6
    theta_b = 0.5
    
    next_theta_a = 0
    next_theta_b = 0
    iter_no = 0
    print '============================================================================================'
    while abs(next_theta_a - theta_a) > 0.00001 and abs(next_theta_b - theta_b) > 0.00001:
        print '-----------------------------------------------------------------------------'
        head_a_sum = 0
        tail_a_sum = 0
        head_b_sum = 0
        tail_b_sum = 0
        if iter_no > 0:
            theta_a = next_theta_a
            theta_b = next_theta_b
        for i in range(n_rounds):
            p_a = (select(n_trials, n_heads[i])) * pow(theta_a, n_heads[i]) * pow(1 - theta_a, n_trials - n_heads[i])
            p_b = (select(n_trials, n_heads[i])) * pow(theta_b, n_heads[i]) * pow(1 - theta_b, n_trials - n_heads[i])
            sum = p_a + p_b
            p_a = p_a / sum
            p_b = p_b / sum
            print 'p_a: %.2f; p_b: %.2f' % (p_a, p_b)
            head_a_sum += n_heads[i] * p_a
            tail_a_sum += (n_trials - n_heads[i]) * p_a
            head_b_sum += n_heads[i] * p_b
            tail_b_sum += (n_trials - n_heads[i]) * p_b
        next_theta_a = head_a_sum / (head_a_sum + tail_a_sum)
        next_theta_b = head_b_sum / (head_b_sum + tail_b_sum)
        print 'Iteration %d: %.4f, %.4f' % (iter_no, theta_a, theta_b)
        print '\t: %.4f, %.4f' % (next_theta_a, next_theta_b)
        iter_no += 1

def simu():
    