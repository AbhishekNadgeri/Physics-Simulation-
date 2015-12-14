#include <iostream>
#include <fstream>

using namespace std;
long double omega;
long double calculate_position(long double position , long double velocity , long double accleration ,long double sampling_time); /* function to calculate position */
long double calculate_velocity(long double velocity , long double accleration , long double sampling_time); /*function to calculate velocity */

int main()
{
     ofstream file;
     file.open("new.dat");
     long double initial_position=0,final_position,initial_velocity,final_velocity,final_time,sampling_time,mass,energy=0,accleration=0,k;
     cout<<"Assuming particle at x=0 at t=0\n"; /* Assuming equation to be y = Asin(wt) */
     cout<<"Enter the value of velocity at t=0 - ";
     cin>>initial_velocity;
     cout<<"Enter the value of Omega  - ";
     cin>>omega;
     cout<<"Enter stop time - ";
     cin>>final_time;
     cout<<"Enter sampling time - ";
     cin>>sampling_time;
     cout<<"enter mass of object -";
     cin>>mass;
     energy = 0.5*(mass*(initial_velocity*initial_velocity)); /* starting there is no inital compression in spring , so all the energy is in the form of kinetic enegry */
     file<<"0"<<" "<<"0"<<" "<<initial_velocity<<" "<<energy<<"\n";
     k= mass*(omega*omega); /* K - is spring constant */
     for(long double t = 0; t< final_time ; t=t+sampling_time)
     {
          final_position = calculate_position(initial_position,initial_velocity,accleration,sampling_time); // function to calaculate final velocity
          final_velocity = calculate_velocity(initial_velocity,accleration,sampling_time);/* here  we find velocity vhich is v(t+/\t/2) */
          accleration = -1*((omega*omega)*final_position);  /* here we find updated accleration a(t+/\t) */
          final_velocity = calculate_velocity(final_velocity,accleration,sampling_time); /* this is the final velocity v(t+/\t) */
          energy = (0.5 *(k* final_position*final_position))+(0.5*mass*(final_velocity*final_velocity)); /* here we find total energy of the system which should be constant */
          file<<t+sampling_time<<" "<<final_position<<" "<<final_velocity<<" "<<energy<<"\n"; /* storing inta a file */
          initial_position = final_position;
          initial_velocity = final_velocity;
     }
     file.close();
     return 0;
}
long double calculate_position(long double position , long double velocity , long double accleration , long double sampling_time) /*Function body to calculate postion of the object */
{
     long double final_position = position + (velocity*sampling_time) + 0.5*(accleration*(sampling_time*sampling_time));
     return final_position;

}

long double calculate_velocity(long double velocity , long double accleration , long double sampling_time) /*Function body to calculate velocity of the object */
{
     long double final_velocity = velocity + 0.5*(accleration*sampling_time);
     return final_velocity;
}


