#### PRINTING RANDOM NUMBERS
import random
from sets import Set
s = Set([])
for x in range(100):
   x = random.randint(1,10111)
   s.add(x)
   #print(random.randint(1,101)),

for i in s:
   print(i),
