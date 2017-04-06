################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Core/Convergence.cpp \
../src/Core/Dictionary.cpp \
../src/Core/GlobalDef.cpp \
../src/Core/InitLBM.cpp \
../src/Core/Parameters.cpp \
../src/Core/Simulation.cpp \
../src/Core/Solution.cpp \
../src/Core/Solver.cpp \
../src/Core/Tau.cpp \
../src/Core/Viscosity.cpp 

OBJS += \
./src/Core/Convergence.o \
./src/Core/Dictionary.o \
./src/Core/GlobalDef.o \
./src/Core/InitLBM.o \
./src/Core/Parameters.o \
./src/Core/Simulation.o \
./src/Core/Solution.o \
./src/Core/Solver.o \
./src/Core/Tau.o \
./src/Core/Viscosity.o 

CPP_DEPS += \
./src/Core/Convergence.d \
./src/Core/Dictionary.d \
./src/Core/GlobalDef.d \
./src/Core/InitLBM.d \
./src/Core/Parameters.d \
./src/Core/Simulation.d \
./src/Core/Solution.d \
./src/Core/Solver.d \
./src/Core/Tau.d \
./src/Core/Viscosity.d 


# Each subdirectory must supply rules for building sources it contributes
src/Core/%.o: ../src/Core/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -o "$@" "$<" -O3 -Wall -c -fmessage-length=0 -lcgns -lhdf5 -lboost_serialization -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)"
	@echo 'Finished building: $<'
	@echo ' '


