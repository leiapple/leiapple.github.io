---
layout: post
title: A start to bash program
categories: Linux
description: bash programming
keywords: linux, bash programming
---

### for loop
```
#/bin/bash
for((i=1;i<=200;i=i+5))
do
echo $i
echo loop.$i.dump
done
```

### if else
```
#/bin/bash
aa=3.14
bb=3.14
one=1
c=`echo "$aa == $bb"|bc`
echo $c
if [[ "$c" -eq 1 ]];then
   echo "ture"
else 
   echo "false"
fi;
```

### if a file exist
```
#this if command is used to judge if a file is empty or not.

if [ -s pidfile ]; then
        echo "hi"
else
        echo "empty"
fi
```

### if string equals

```
flag_0=DISLOCATION

if [ "$flag_0" = "DISLOCATION" ]; then
echo "equal 000"
else
echo "not equal 000"
fi
```