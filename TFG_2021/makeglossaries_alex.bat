set path=%path%;"Z:\Programs\MiKTeX\miktex\bin\x64"

makeindex memoria_TFM.idx

makeindex -s "memoria_TFM.ist" -t "memoria_TFM.alg" -o "memoria_TFM.acr" "memoria_TFM.acn"

makeindex -s "memoria_TFM.ist" -t "memoria_TFM.glg" -o "memoria_TFM.sym" "memoria_TFM.sbl"

makeindex -s "memoria_TFM.ist" -t "memoria_TFM.glg" -o "memoria_TFM.gls" "memoria_TFM.glo"
