# python verify_q11.py -k 3 -l 100 -i input/k_3_n_100.npz -s input/k_3_n_100.fna -o output/k_3_n_100.fna
# python verify_q11.py -k 2 -l 100 -i input/random_k2.npz -o output/random_k2.fna
# python verify_q11.py -k 5 -l 10000 -i input/long.npz -s input/long.fna -o output/long.fna
for k in 1 2 3 4 5; do
  for l in 50 100 200 500 1000 10000; do
    echo "k: ${k}, l: ${l}"
    python verify_q11.py -k $k -l $l -i "input/random_prob_k${k}.npz" -o "output/random_prob_k_${k}_n_${l}.fna"
  done
done
