################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Algorithm/LowOrder/Boundaries.cpp \
../src/Algorithm/LowOrder/CollideLowOrder.cpp \
../src/Algorithm/LowOrder/D2Q9.cpp \
../src/Algorithm/LowOrder/D2Q9ColourFluid.cpp \
../src/Algorithm/LowOrder/D2Q9_TwoPhases.cpp \
../src/Algorithm/LowOrder/Gradients.cpp \
../src/Algorithm/LowOrder/GradientsDEF.cpp \
../src/Algorithm/LowOrder/StreamLowOrder.cpp 

OBJS += \
./src/Algorithm/LowOrder/Boundaries.o \
./src/Algorithm/LowOrder/CollideLowOrder.o \
./src/Algorithm/LowOrder/D2Q9.o \
./src/Algorithm/LowOrder/D2Q9ColourFluid.o \
./src/Algorithm/LowOrder/D2Q9_TwoPhases.o \
./src/Algorithm/LowOrder/Gradients.o \
./src/Algorithm/LowOrder/GradientsDEF.o \
./src/Algorithm/LowOrder/StreamLowOrder.o 

CPP_DEPS += \
./src/Algorithm/LowOrder/Boundaries.d \
./src/Algorithm/LowOrder/CollideLowOrder.d \
./src/Algorithm/LowOrder/D2Q9.d \
./src/Algorithm/LowOrder/D2Q9ColourFluid.d \
./src/Algorithm/LowOrder/D2Q9_TwoPhases.d \
./src/Algorithm/LowOrder/Gradients.d \
./src/Algorithm/LowOrder/GradientsDEF.d \
./src/Algorithm/LowOrder/StreamLowOrder.d 


# Each subdirectory must supply rules for building sources it contributes
src/Algorithm/LowOrder/%.o: ../src/Algorithm/LowOrder/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -O0 -g3 -Wall -c -fmessage-length=0 -lcgns -lhdf5 -lboost_serialization -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


