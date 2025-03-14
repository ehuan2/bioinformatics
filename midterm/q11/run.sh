# python q11.py --k 2 --L 100 --p input/random_k2.npz
# python q11.py --k 5 --L 10000 --p input/long.npz
for k in 1 2 3 4 5; do
  for n in 50 100 200 500 1000 10000; do
    echo "Running k: ${k}, n: ${n}"
    echo ">seq0" > "output/random_prob_k_${k}_n_${n}.fna"
    python q11.py --L $n --k $k --p "input/random_prob_k${k}.npz" >> "output/random_prob_k_${k}_n_${n}.fna"
  done
done
