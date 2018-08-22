calc () {
    bc -l <<< "$@"
}

export NUM=$(grep "Iteration" out | grep "of" | awk '{print $6}' | wc -l)
export SUM=$(grep "Iteration" out | grep "of" | awk '{print $6}' | awk '{s+=$1} END {print s}')

export SUMMA=$(sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' <<<"$SUM")

echo $NUM, $SUMMA
export MEAN=$(calc $SUMMA/$NUM)
echo "Time for iteration: " $MEAN

export NUM=$(grep "Converged after" debug_info.txt | awk '{print $3}' | wc -l)
export SUM=$(grep "Converged after" debug_info.txt | awk '{print $3}' | awk '{s+=$1} END {print s}')
echo $NUM, $SUM
export MEAN=$(calc $SUM/$NUM)
echo "Converged after: " $MEAN


