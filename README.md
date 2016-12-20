# ConjugateGradientsMethod
## Построение графиков
Построение графиков производится с помощью системы Matlab
## Сборка
module add   impi/5.0.1
module add   intel/15.0.090
module add   cuda/6.5.14
nvcc -ccbin mpicxx ConjugateGradientsClass.cpp ConjugateGradientsMethod.cu  GridClass.cpp GridCLassCUDA.cu InputFunctions.cpp Main.cu  XGetopt.cpp -o CGM
sbatch -p gputest -n 4 impi ./CGM

## Запуск
### Аргументы
-x int -- нижний левый угол сетки ось X

-y int  -- нижний левый угол сетки ось Y

-a int -- верхний правый угол сетки ось X

-b int -- верхний правый угол сетки ось Y

-r int -- число строк в сетке

-c int -- число столбцов в сетке

-f string -- выходной файл

Функции F и phi задаются в файле InputFunctions.сu и далее используются как callback-функции
На данный момент в Main.cu заданы тестовые значения, поэтому аргументы командной строки не требуютс

### На СК Ломоносов
sbatch -p gputest -n 4 impi ./CGM
