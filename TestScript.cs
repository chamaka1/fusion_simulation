using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.EventSystems;
using TMPro;



public class TestScript : MonoBehaviour
{
    
    public Rigidbody rb;
    
    public Vector3 sphere_position;
    public double mu;
    public const int particle_no = 1;
    public double ds; // spacing between each particle
    public double q; // charge
    public double m; // mass
    public double c; // speed of light
    public double I_coils; // current in coils
    public double N;
    public double I_plasma;
    public double g;
    public const double a = 1.5; // radius of each coil
    public const double b = 0.8; // radius of central region
    public double t_final; // duration of sim
    public double dt; // step size


    public double[] xp = new double[particle_no];
    public double[] yp = new double[particle_no];
    public double[] zp = new double[particle_no];
    public double[] vxp = new double[particle_no];
    public double[] vyp = new double[particle_no];
    public double[] vzp = new double[particle_no];
    public double[] gammap = new double[particle_no];

    // initial condition
    public double[] initial_position = new double[3];
    public double[] initial_speed = new double[3];


    
    // Start is called before the first frame update
    void Start()
    {
        double[] initial_position = {2,-1,0} ;
        double[] initial_speed = {-0.1,-0.15,0} ;
        // sim conditions
        t_final = 1.2f * Math.Pow(10, -5);
        dt = 1f * Math.Pow(10, -10); 
        //particle_no = 1 ;
        ds = 3/particle_no; // spacing diameter of whole coil
        
                

        for (int i = 0; i < particle_no ; i++)
        {
            xp[i] = initial_position[0];
            yp[i] = initial_position[1];
            zp[i] = ds * particle_no - ds*(i+1);
            vxp[i] = initial_speed[0];
            vyp[i] = initial_speed[1];
            vzp[i] = initial_speed[2];

        }
        
        // environment variables
        mu = 4 * Math.PI * Math.Pow(10, -7) ;
        q = 1.60217 * Math.Pow(10,-19);
        m = 1.67262 * Math.Pow(10, -27);
        c = 3 * Math.Pow(10, 8);
        g = 9.81;


        rb = GetComponent<Rigidbody>();
        
        sphere_position = new Vector3((float) xp[0], (float)yp[0],(float) zp[0]);
        GetComponent<Rigidbody>().position = sphere_position;
        //Rigidbody.MovePosition(position);
        
       
        Debug.Log(sphere_position);
    }

    // fixed update must be set to be same as dt
    void FixedUpdate()
    {
        for (int i = 0; i < particle_no ; i++) // loop for each particle
        {

            double x = xp[i];
            double y = yp[i];
            double z = zp[i];
            double vx = vxp[i];
            double vy = vyp[i];
            double vz = vzp[i];
            double gamma = gammap[i];

            double theta = Math.Atan2(y, x);
            //distance = sqrt( z^2 + (x-(a+b)*cos(phi))^2 + (y-(a+b)*sin(phi))^2 );
            double[] input = new double[]{x, y, z};
            double[] B = toroid(input);

            double Bx = B[0];
            double By = B[1];
            double Bz = B[2];

            double ax = q/(gamma*m) * (vy*Bz-vz*By);
            double ay = q/(gamma*m) * (vz*Bx-vx*Bz); 
            double az = q/(gamma*m) * (vx*By-vz*Bx) -g;
            Debug.Log(ax);
            vx = vx + ax*dt;
            vy = vy + ay*dt;
            vz = vz + az*dt;
            //Debug.Log(x);
            //Debug.Log(vx);  

            //sim error
            vxp[i] = vx * Math.Sqrt(Math.Pow(initial_speed[0], 2) + Math.Pow(initial_speed[1], 2) + Math.Pow(initial_speed[2], 2)) / (Math.Sqrt(Math.Pow(vx, 2) + Math.Pow(vy, 2) + Math.Pow(vz, 2))) / Math.Pow(c, 2);
            vyp[i] = vy * Math.Sqrt(Math.Pow(initial_speed[0], 2) + Math.Pow(initial_speed[1], 2) + Math.Pow(initial_speed[2], 2)) / (Math.Sqrt(Math.Pow(vx, 2) + Math.Pow(vy, 2) + Math.Pow(vz, 2))) / Math.Pow(c, 2);
            vzp[i] = vz * Math.Sqrt(Math.Pow(initial_speed[0], 2) + Math.Pow(initial_speed[1], 2) + Math.Pow(initial_speed[2], 2)) / (Math.Sqrt(Math.Pow(vx, 2) + Math.Pow(vy, 2) + Math.Pow(vz, 2))) / Math.Pow(c, 2);


            xp[i] = x + vxp[i]*dt;
            yp[i] = y + vyp[i]*dt;
            zp[i] = z + vzp[i]*dt;
         
        }

        sphere_position = new Vector3((float) xp[0], (float) yp[0], (float) zp[0]);
        GetComponent<Rigidbody>().position = sphere_position;
        //Debug.Log(sphere_position);



    }

    public double[] toroid(double[] input) {
        double[] B = new double[]{0, 0, 0};

        double r = Math.Sqrt(Math.Pow((double)input[0], 2) + Math.Pow((double)input[1], 2));
        double B_r = -2f * (1 / r) * input[2] * Math.Pow(r,2);
        double B_z = -1f / r *((4f * r - 4f * Math.Pow(r , 3)) - (2 * r * Math.Pow(input[2] , 2)));
        double B_theta = 0f;

        double theta = Math.Atan2(input[1], input[0]);

        B[0] = Math.Cos(theta) * B_r - Math.Sin(theta) * B_theta;
        B[1] = Math.Sin(theta) * B_r + Math.Cos(theta) - B_theta;
        B[2] = B_z;

        return B;
    }   
}

