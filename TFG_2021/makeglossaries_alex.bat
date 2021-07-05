set path=%path%;"Z:\Programs\MiKTeX\miktex\bin\x64"

REM makeindex memoria_TFG.idx

makeindex -s "memoria_TFG.ist" -t "memoria_TFG.alg" -o "memoria_TFG.acr" "memoria_TFG.acn"

REM makeindex -s "memoria_TFG.ist" -t "memoria_TFG.glg" -o "memoria_TFG.sym" "memoria_TFG.sbl"

makeindex -s "memoria_TFG.ist" -t "memoria_TFG.glg" -o "memoria_TFG.gls" "memoria_TFG.glo"

pause