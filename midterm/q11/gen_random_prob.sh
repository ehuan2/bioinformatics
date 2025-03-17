for k in {1..5}; do
    python gen_random_prob.py -k $k --output input/random_prob_k$k
done
