# Technical tricks

## Loading software

After logging in, and also after having started a screen below, do the
following:

`source IN-BIOS9000-software-load`



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
