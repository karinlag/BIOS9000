# Technical tricks

## Loading software

After logging in, you will need software.

See the list here:

https://uio-in-biosx000.readthedocs.io/en/latest/software.html

You will have to load the individual modules per software to run.

## Getting a compute session

Quite frequently we will have to use slurm to run things. 

The command to get an interactive session is the following:

'salloc --ntasks=4 --mem-per-cpu=4G --time=03:00:00 --qos=devel --account=ec34'

This will give you 4 cpus for 3 hours, with 20GB of memory.

## Using screen

`screen` is a small program that will allow programs to keep running even if
you close your terminal or lose internet access.

### How to start

1. Type in `screen`. You will get a new prompt, but that's it
2. Start whatever program that yo want to have running.
3. Depress `Control-a`, followed by `d`. Your screen will now say something
about a screen being `detached`.

You can now end your ssh session, and your program will keep running.

### How to see and resume your sessions

You can see a list of screens (you can have multiple), type in `screen -ls`.
You will see something along the lines of `23921.pts-123.login-0-1`. The number
before `pts` is what I refer to below.
To connect to a specific screen, type in `screen -DR <number>`.

For more, [check out this link](https://www.tecmint.com/screen-command-examples-to-manage-linux-terminals/)

## Using nodes

For doing the actual assemblies, we will ask you to log into some interactive
compute nodes to run stuff. This is done by using a command called `ssh`, like
this:

`ssh int-X`

where X is a number. 

You will then be asked for the authentication number again, and your password.




