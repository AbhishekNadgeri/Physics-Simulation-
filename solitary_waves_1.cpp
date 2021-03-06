#include <iostream>		//This is a code to try to see how solitay wave propogate through force and velocity dependencies 
#include<math.h>		//After implementing code a 'new.dat' file will be created and can be plotted using gnu plot or matlab 
#include<fstream>
#include<stdio.h>

using namespace std;
int n;
class sphere{
long double  R,V,E,U,velocity,accleration,mass;   // R - Radius , V - poisson's ratio , E - Young's Modulus , U - current position of particle
public:
void set_values(long double  radius,long double  poissons_ratio,long double  youngs_modulus,long double  position,long double  m);
long double  get_radius();
long double  get_position();
long double  get_poissons_ratio();
long double  get_youngs_modulus();
long double  get_velocity();
long double  get_accleration();
long double  get_mass();
long double  set_position(long double  pos);
long double  set_velocity(long double  vel);
long double  set_accleration(long double  accl);
long double  set_mass(long double  m);
};
void sphere::set_values(long double  radius,long double  poissons_ratio,long double  youngs_modulus,long double  position,long double  m) //Constructor for setting the values
{
    R=radius;
    V=poissons_ratio;
    E=youngs_modulus;
    U=position;
    velocity=0;
    accleration=0;
    mass=m;
}

//returns radius of sphere
long double  sphere::get_radius()
{
    return R;
}

//returns position of particle
long double  sphere::get_position()
{
    return U;
}

//returns poissons ratio
long double  sphere::get_poissons_ratio()
{
    return V;
}

//return youngs modulus
long double  sphere::get_youngs_modulus()
{
    return E;
}

//returns velocity of particle
long double  sphere::get_velocity()
{
    return velocity;
}

//returns particle accleration
long double  sphere::get_accleration()
{
    return accleration;
}

//returns particle mass
long double  sphere::get_mass()
{
    return mass;
}

//sets position of particle
long double  sphere::set_position(long double  pos)
{
    U=pos;
}

//sets velocity of particle
long double  sphere::set_velocity(long double  vel)
{
    velocity=vel;
}

//sets accleration of particle
long double  sphere::set_accleration(long double  accl)
{
    accleration=accl;
}

//sets mass of particle
long double  sphere::set_mass(long double  m)
{
    mass=m;
}

//to check if the two particle are in contact
bool in_contact (sphere s1 , sphere s2)
{
    long double position = s1.get_position() - s2.get_position();
    if(position<0)
        position = position * -1;
    if(position <= (s1.get_radius() + s2.get_radius()) )
        return true;
    return false;
}

//calculate value of D for the particle
long double  calculate_D (sphere s1, sphere s2)
{
long double  d1,d2,D;
d1 = ( 1 - ( s1.get_poissons_ratio() * s1.get_poissons_ratio() ) ) / s1.get_youngs_modulus();
d2 = ( 1 - ( s2.get_poissons_ratio() * s2.get_poissons_ratio() ) ) / s2.get_youngs_modulus();
D = 0.75 * ( d1 + d2 );
return D;
}

//calculates value of A for the particle
long double  calculate_A(sphere s1,sphere s2)
{

    long double  a = sqrt(((s1.get_radius()*s2.get_radius()) / (s1.get_radius()+s2.get_radius())) );
    long double  d = calculate_D(s1,s2);
    long double  A = (0.4*a) / d;
    A = 2.5 * A;
    return A;
}

//calculates position of the particle
long double  calculate_position(sphere *s1,int i,long double  sampling_time)
{
    long double  final_position = s1[i].get_position() +  (s1[i].get_velocity() * sampling_time) + 0.5*(s1[i].get_accleration()*sampling_time*sampling_time); // s(t+/\t) = s(t) + (v * /\t) + (0.5 * a * /\t * /\t)
    return final_position;
}

//calculates velocity of particle
long double  calculate_velocity(sphere s1,long double  sampling_time)                             //using velocity verlet
{

    long double  final_velocity = s1.get_velocity() + 0.5*(s1.get_accleration() * sampling_time); // v = u + 0.5 * a * /\t
    return final_velocity;
}

//calculates accleration of particle
long double  calculate_accleration(sphere *s1, int i)
{
    long double  accleration,temp;                                  // a = A(i-1,i) * [ (r(i-1) + r(i) + u(i-1) - u(i)) ^ 1.5 ] - A(i,i+1) * [ (r(i) + r(i+1) + u(i) - u(i-1) ]
    if(i==0)
    {
        if(in_contact(s1[0],s1[1])==true)
        {
           temp = s1[0].get_position() - s1[1].get_position();
            if(temp<0)
                temp=temp*-1;
            temp = s1[0].get_radius() + s1[1].get_radius() - temp;
            accleration = -1 * calculate_A(s1[0],s1[1]) * pow(temp, 1.5);
            accleration = accleration / s1[0].get_mass();

        }
        else
        accleration = 0;
    }
    else if(i==n-1)
    {
        if(in_contact(s1[n-2],s1[n-1])==true)
        {
            temp = s1[n-2].get_position() - s1[n-1].get_position();
            if( temp < 0 )
                temp=temp*-1;
            temp = s1[n-2].get_radius() + s1[n-1].get_radius() - temp;
            accleration = calculate_A(s1[n-2],s1[n-1]) * pow(temp , 1.5);
            accleration = accleration / s1[n-1].get_mass();
        }
        else
        accleration=0;
    }
    else
    {

        long double  a1,a2;
        if(in_contact(s1[i-1],s1[i])==true)
        {
            temp = s1[i-1].get_position() - s1[i].get_position();
            if( temp < 0 )
                temp = temp * -1;
            temp = s1[i-1].get_radius() + s1[i].get_radius() - temp;
            a1 = calculate_A(s1[i-1],s1[i]) * pow(temp , 1.5);
            a1 = a1 / s1[i].get_mass();
        }
        else
            a1 = 0;
        if(in_contact(s1[i],s1[i+1])==true)
        {
            temp = s1[i].get_position() - s1[i+1].get_position();
            if( temp < 0)
                temp = temp * -1;
            temp = s1[i].get_radius() + s1[i+1].get_radius() - temp;
            a2 = calculate_A(s1[i],s1[i+1]) * pow(temp , 1.5);
            a2 = a2 / s1[i].get_mass();
        }
        else
            a2 = 0;
        accleration = a1 - a2;

    }



    return accleration;

}


int main()
{

    long double  radius,poisson,young,mass,sampling_time,simulation_time,striker,density;
    int p;
    cout<<"Enter the number of particels in the system - ";
    cin>>n;
    sphere s[n];
    cout<<"Assuming that the particles are of same material and initally all particles are touching each other . Enter the following details \n";
    cout<<"Enter the radius of particles - ";
    cin>>radius;
    cout<<"Enter Poisson's ratio of the material sphere - ";
    cin>>poisson;
    cout<<"Enter Young's Modulus of the material sphere - ";
    cin>>young;
    cout<<"Enter density of particles - ";
    cin>>density;
    cout<<"Enter Sampling time - ";
    cin>>sampling_time;
    cout<<"Enter Simulation time (Assuming is starts from t=0) - ";
    cin>>simulation_time;
    mass =((4 * 3.14159265359 * radius * radius * radius)/3)*density;


    for(int i=0;i<n;i++)
    {
        s[i].set_values(radius , poisson , young , 2*i*radius , mass); // setting all the respective values of particle
    }

    cout<<"Enter the velocity of the striker particle - ";
    cin>>striker;
    s[0].set_velocity(striker);
    cout<<"Enter the nth particle for which you want the data - ";
    cin>>p;


    ofstream file;                  //opening file to write
    file.open("new.dat");



    for(long double  i=0;i<=simulation_time;i=i+sampling_time)  // loop for calculating values of all particle for all time period
    {
        file<<i<<" ";
        if(i==0)
        file<<"0 "<<s[p-1].get_velocity()<<" "<<s[p-1].get_position();
        else
        {
            for(int j=0;j<n;j++)                               // loop for calculating values of a specific particle at time i
            {
                long double  final_position = calculate_position(s,j, sampling_time);       // calculating position r(t+/\t)
                long double  final_velocity = calculate_velocity(s[j] , sampling_time);     // calculating velocity v(t+/\t/2)
                long double  accleration = calculate_accleration(s,j);                      // calculating accleration a(t+/\t)
                s[j].set_velocity(final_velocity);                                          // sets velocity of particle v(t+/\t/2)
                s[j].set_accleration(accleration);                                          // sets accleration of particle
                final_velocity = calculate_velocity(s[j] , sampling_time);                  // calculating velocity v(t+/\t)
                s[j].set_velocity(final_velocity);                                          // sets velocity of particle v(t+/\t)
                s[j].set_position(final_position);                                          // sets position of particle r(t+/\t)
                long double  force = s[j].get_mass() * s[j].get_accleration();              // calculates force acting on paticle
                if((j+1)==p)
                    file<<force<<" "<<s[j].get_velocity()<<" "<<s[j].get_position();        //writing values into the file
            }
        }

        file<<"\n";

    }
    file.close();
    return 0;
}



