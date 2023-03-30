###############################################################
# Introduction to R                                           #
# By: ArrietaLab - University of Calgary                      #
# Author: Erik van Tilburg Bernardes                          # 
# Dates: 3-5 April 2023                                       #
# Location: Oswaldo-Cruz Institute - Rio de Janeiro (Brazil)  #
###############################################################

#Create object "x"
x <- 5
x #you just assigned value "5" to an object called "x"


#create additional objects
#objects can be made up of alphanumeric characters but can't start with a number. 
#Object names are case sensitive.

y <- 2
students <- 20
A <- "Fiocruz"
y
#print the values for the different objects created


#Operations in R

#There are numerous operators in R - arithmetic, logical tests, etc.
# you have created objects x and y before
# now create an object that is the sum of them
z <- x + y
z

#or you can just use them directly for different operations 
x * y # value if you multiply them

# let's use some logic now
x > y

#what about x and z?


#Functions

# let's run function "ls()"
ls() #what this function does is to list all the objects and functions created in your environment


#let's explore the "help()" function
# on it's own it displays documentation for itself...try it...
help() 

#you can also add the name of another function to explore what it does
#let's get more details on the "ls()" function we just used
help(ls)  # same as '?ls'

#let's explore other built-in functions
mean(x, y, z)
min(x, y, z)
print(x + y + z)
sum(x, y, z)


#let's create our own function
MySequence <- seq(1, 100)
MySequence #print your new created sequence

MyMean <- mean(MySequence)
MyMean

SumMySequence <- sum(MySequence)
SumMySequence
print(sum(MySequence))


#vector objects
# a vector is the simplest type of object in R.
# it is an object that contains one or more values of the same type. Some common types of vector are:

vec.a <- c(x, y, z)
vec.a

#so now creating a new vector
vec.b <- c("x", "y", "z")
vec.b #is it the same? what is different?

#Matrix objects
#A matrix is a two-dimensional numerical array
#let's create a matrix
mat.a <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
mat.a

#let's name the rows and columns of the new matrix
rownames(mat.a) <- c("row1", "row2")
mat.a
colnames(mat.a) <- c("col1", "col2", "col3")
mat.a

#data.frame objects
#its an object, similar to matrix, but its columns can be a mixture of different types of data. 
# The data.frame is perhaps the most frequently used type of object for biodiversity analysis since it can contain many different types of data.

df.a <- as.data.frame(mat.a)
df.a$col4 <- c("foo", "bar")
dim(df.a)
## [1] 2 4
df.a

#List objects
#A list is a collection of different types of objects. 
# make a list combining the vector, matrix, and data.frame objects we just created
list.a <- list(va = vec.a, vb = vec.b, ma = mat.a, dfa = df.a)
names(list.a)
list.a



#Tidyverse

#let's install this new package
library(tidyverse) #note that tidyverse comes wiht several other packages attached to it
help(tidyverse)


#import global variation in atopy and asthma file
getwd()
atopyasthmavariation <- read.table("WorldwideVariationAtopy.tsv", header=TRUE, sep="\t", dec=".", strip.white=TRUE)
atopyasthmavariation
# read.table is a function to read a data table in 'tsv' format that we have in our workspace
# part of the function is to tell if the table has a header (=TRUE)
# sep indicated field separator between characters
# dec specifies the decimal points used in the file (it might be using ",")
# strip.white required as we specified "sep" that will account for blank spaces between characters


# Let's explore the new imported data.frame
colnames(atopyasthmavariation)
dim(atopyasthmavariation)


# too much information,  what can we do?
#let's filter and sort the data:

AAvariationEurope <- filter(atopyasthmavariation, 
                            Continent == "Europe")
AAvariationEurope
dim(AAvariationEurope)
# you used the equal to sign "==" to select only the data you wanted

# what is you want to do the oposite? instead of selecting just one, you want to exclude a continent
AAvariationNoOceania <- filter(atopyasthmavariation, 
                               Continent != "Oceania")
AAvariationNoOceania
dim(AAvariationNoOceania)
# you used the differ to sign "!=" to select all  data but the one you dont want

#Now we want to see the countries with the highest atopy incidence
HighAtopy <- filter(atopyasthmavariation, 
                    Atopy > 25) %>% 
  arrange(Country)
HighAtopy
# you now selected only the countries with Atopy values higher than 25% of the sample population
# here we introduce a new language: the pipe operator in the form %>% aims to combine various functions without the need to assign the result to a new object. The pipe operator passes the output of a function applied to the first argument of the next function. 
# using the pipe operator, we applied the filter and the arrange function to the same object HighAtopy

#we can also group informations 
AApercontinent <- atopyasthmavariation %>% 
  group_by(Continent) %>% 
  summarise(mx = max(Asthma_Ever, na.rm = TRUE), 
            min = min(Asthma_Ever, na.rm = TRUE),
            median = median(Asthma_Ever))

AApercontinent
# now we created an object AApercontinent, in which we grouped the Asthma_Ever incidence per continent
# in addition, we summarized the data for max, minimum, and median values.


#Now lets see the data!

#lets start with some built-in tools
hist(atopyasthmavariation$Asthma_Ever)

#what about atopy?
hist(atopyasthmavariation$Atopy)

#now let's do more...
#first install ggplot2
library(ggplot2)
#check if it was properly installed
#also inspect ggplot2 to get more information about the package

plot1 <- ggplot(atopyasthmavariation, aes(Asthma_Ever, Atopy, colour = Continent)) + 
  geom_point()
plot1
# here we are just plotting our data
# asthma_ever vs. atopy, and color coded per continent
# but as is we can't make much sense of it


plot2 <- ggplot(atopyasthmavariation, aes(Asthma_Ever, Continent, colour = Continent)) + 
  geom_point()
plot2
# here we are now looking at asthma incidence per continent

plot3 <- ggplot(atopyasthmavariation, aes(Asthma_Ever, Continent, colour = Continent)) + 
  geom_boxplot()
plot3
# here we have the same graph...but instead of geom_points, now we are using geom_boxplot
# asthma_ever vs. atopy, and color coded per continent

# but we mentioned with ggplot we can edit the graphs so they look however we want, we just keeping adding layers to it...so let's explore
plot4 <- ggplot(atopyasthmavariation, aes(Asthma_Ever, Continent, fill = Continent, colour = Continent)) + #adding fill, now we also colored the boxplots instead of just outlines
  geom_boxplot(alpha=1) + #keeping the assigned boxplots for the graph
  geom_jitter(width=0.2) + #but we also want to see the individual points, not just the bars
  xlab("Asthma Incidence (%)") + #now we are renaming the x axis
  scale_color_manual(name="Continent", values=c("black", "black", "black", "black", "black", "black", "black")) 
plot4

# now let's do the same for atopy
plot5 <- ggplot(atopyasthmavariation, aes(Atopy, Continent, fill = Continent, colour = Continent)) + 
  geom_boxplot(alpha=1) +
  geom_jitter(width=0.2) + 
  xlab("Atopy Incidence (%)") + #here we need to manually change the axis name
  scale_color_manual(name="Continent", values=c("black", "black", "black", "black", "black", "black", "black")) 
plot5

# and another nice thing we can do is to combine graphs to generate a report
# we need to install a new package though
library(gridExtra)
grid.arrange(nrow = 2,  plot4, plot5)
