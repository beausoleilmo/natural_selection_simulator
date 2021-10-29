# Natural selection simulator

In this repository, I explore how to program natural selection in R. 

- Bacteria are generated in an environment with food units
- Bacteria are moving randomly in the environment 
- When a bacteria is "in range" of a food unit, it can eat it: it gains fitness units 
- After X amount of food units = survives, and after Y amount of food unit, it reproduces 
- When bacteria reproduce, the have a progeny that have the parental traits with some mutation (random fluctution in the phenotype). 


The script can generate plots (from base R) and a gif (using [ImageMagick](https://imagemagick.org/index.php)). 


<!-- ![Example of a simulation run by the program](gif/ns.film.gif)-->
<img src="gif/ns.film.gif" width="400" height="400" />

The bacteria are in yellow, the food in red. 


# Things to modify or add 

- Cost of having higher speed? 
  - i.e., costs more to go faster or the energy spent is proportional to the speed
- Add a fitness function
- try resetting the position of bacteria and of food after search time is done 


# Including more traits 

# Inspirations 
- https://www.youtube.com/watch?v=0ZGbIKd0XrM&ab_channel=Primer 
- http://www.netlogoweb.org/launch#http://www.netlogoweb.org/assets/modelslib/Curricular%20Models/GenEvo/GenEvo%203%20Genetic%20Drift%20and%20Natural%20Selection.nlogo
