#!/bin/bash
if [ -d force_x ];
then
  rm -f force_x/*
  echo "rm -f force_x/*"
fi

if [ -d force_y ];
then
  rm -f force_y/*
  echo "rm -f force_y/*"
fi
