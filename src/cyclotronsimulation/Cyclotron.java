package cyclotronsimulation;

import java.awt.Color;
import java.awt.Graphics;
import java.util.ArrayList;

/**
 *
 * @author Johannes Bechberger
 */
public class Cyclotron extends javax.swing.JPanel {

   private double radius_m; //in m
   private double px_per_m;
   private double s = 0; //Distance between the to half circles in m
   private double magnetic_induction;
   private double current;
   private double objects_per_second;
   private ArrayList<FlyingObject> flying_objects = new ArrayList<FlyingObject>();
   private int object_radius = 10; //in px
   private double object_mass; //in proton masses
   private long last_object_added_time = 0; //in ms
   private long last_time = -1; //in ms
   private double max_v = 0;
   private double max_n = 0;
   private double max_r = 0;
   private double time_fac = 0;

   /**
    * @return the objects_per_second
    */
   public double getObjects_per_second() {
      return objects_per_second;
   }

   /**
    * @param objects_per_second the objects_per_second to set
    */
   public void setObjects_per_second(double objects_per_second) {
      this.objects_per_second = objects_per_second;
   }

   /**
    * @return the flying_objects
    */
   public ArrayList<FlyingObject> getFlying_objects() {
      return flying_objects;
   }

   /**
    * @param flying_objects the flying_objects to set
    */
   public void setFlying_objects(ArrayList<FlyingObject> flying_objects) {
      this.flying_objects = flying_objects;
   }

   /**
    * @return the object_radius
    */
   public int getObject_radius() {
      return object_radius;
   }

   /**
    * @param object_radius the object_radius to set
    */
   public void setObject_radius(int object_radius) {
      this.object_radius = object_radius;
   }

   /**
    * @return the object_mass
    */
   public double getObject_mass() {
      return object_mass;
   }

   /**
    * @param object_mass the object_mass to set
    */
   public void setObject_mass(double object_mass) {
      this.object_mass = object_mass;
      for (FlyingObject flyingObject : flying_objects) {
         flyingObject.setMass(object_mass);
      }
   }

   /**
    * @return the last_object_added_time
    */
   public long getLast_object_added_time() {
      return last_object_added_time;
   }

   /**
    * @param last_object_added_time the last_object_added_time to set
    */
   public void setLast_object_added_time(long last_object_added_time) {
      this.last_object_added_time = last_object_added_time;
   }

   /**
    * @return the last_time
    */
   public long getLast_time() {
      return last_time;
   }

   /**
    * @param last_time the last_time to set
    */
   public void setLast_time(long last_time) {
      this.last_time = last_time;
   }

   /**
    * @return the time_fac
    */
   public double getTime_fac() {
      return time_fac;
   }

   /**
    * @param time_fac the time_fac to set
    */
   public void setTime_fac(double time_fac) {
      this.time_fac = time_fac;
   }

   public class FlyingObject {

      private double time = 0; //in s
      private double velocity = 0; //in m/s
      private double last_velocity = 0; //in m/s
      private double mass; //in kg
      private double acceleration; //in m/s²
      private double trace_radius = 0;
      private double last_trace_radius = 0;
      private double last_last_trace_radius = 0;
      private double trace_time;
      //(0m|0m) is the middle of the cyclotron
      private double x = 0; //in m
      private double y = 0; //in m
      private double last_x = -1;
      private double last_y = -1;
      private double last_in_y = -1;
      private double n = 4;
      private double in_middlearea_time = 0;
      private double in_half_circle_time = 0;
      private boolean was_in_middlearea = true;
      public static final double Q = 1.60218E-19; //in C
      public static final double PROTON_MASS = 1.67E-27; //in kg

      /**
       *
       * @param object_radius in px
       * @param mass in proton masses
       */
      public FlyingObject(double mass) {
         this.mass = mass * PROTON_MASS;
      }

      public void setCoordinates(double x, double y) {
         this.x = x;
         this.y = y;
      }

      public boolean isInMiddleSpace() {
         return x <= s / 2 && x >= -s / 2 && y <= radius_m && y >= -radius_m;
      }

      public boolean isOutOfZyklotronCoreArea() {
         return isOutOfHalfCircles() && !isInMiddleSpace();
      }

      public boolean isOutOfHalfCircles() {
         return Math.sqrt(Math.pow(Math.abs(x) - (s / 2), 2) + Math.pow(y, 2)) > radius_m && !isInMiddleSpace();
      }

      public boolean isOutOfZyklotron() {
         return Math.abs(x * px_per_m) > getWidth() / 2 || Math.abs(y * px_per_m) > getHeight() / 2;
      }

      public void updateCoordinates(double add_time) {
         time += add_time;
         double _x = x;
         double _y = y;
         if (isOutOfZyklotronCoreArea()) {
            x = x + (x - last_x);
            y = y + (y - last_y);
         } else {
            if (isInMiddleSpace()) {
               acceleration = (n % 2 == 0 ? -1 : 1) * (n == 0 ? 0.5 : 1) * current * Q / (mass * s);
               double timespan = time - in_middlearea_time;
               velocity = (n % 2 == 0 ? -1 : 1) * (Math.abs(acceleration * timespan) + Math.abs(last_velocity));
               double x_span = 0.5 * acceleration * Math.pow(timespan, 2) + last_velocity * timespan;
               x = (n == 0 ? 0 : (x_span > 0 ? -s / 2.0 : s / 2.0)) + x_span;
               was_in_middlearea = true;
               in_half_circle_time = time;
               last_in_y = y;
            } else {
               in_middlearea_time = time;
               if (was_in_middlearea) {
                  n++;
                  System.out.println("velocity = " + velocity);
                  was_in_middlearea = false;
                  last_velocity = velocity;
                  last_trace_radius = last_last_trace_radius;
               }
               trace_radius = Math.abs(mass * velocity / (Q * magnetic_induction));
               last_last_trace_radius = trace_radius;
               trace_time = 2 * Math.PI * (mass / Q) / magnetic_induction;
               double timespan = time - in_half_circle_time;
               x = (n % 2 == 0 ? -1 : 1) * trace_radius * Math.sin(2 * Math.PI / trace_time * timespan) + (n % 2 == 0 ? -s / 2.0 : s / 2.0);
               y = (n % 2 == 0 ? -1 : 1) * (trace_radius * Math.cos(2 * Math.PI / trace_time * timespan) + (Math.abs(last_in_y) - trace_radius));
            }
         }
         max_n = max_n < n ? n : max_n;
         max_v = max_v < Math.abs(velocity) ? Math.abs(velocity) : max_v;
         max_r = max_r < Math.abs(trace_radius) ? Math.abs(trace_radius) : max_r;
         last_x = _x;
         last_y = _y;
      }

      public void setMass(double mass) {
         this.mass = mass * PROTON_MASS;
      }

      /**
       * @return the velocity
       */
      public double getVelocity() {
         return velocity;
      }

      /**
       * @param velocity the velocity to set
       */
      public void setVelocity(double velocity) {
         this.velocity = velocity;
      }

      /**
       * @return the mass
       */
      public double getMass() {
         return mass;
      }

      /**
       * @return the acceleration
       */
      public double getAcceleration() {
         return acceleration;
      }

      /**
       * @param acceleration the acceleration to set
       */
      public void setAcceleration(double acceleration) {
         this.acceleration = acceleration;
      }

      /**
       * @return the x
       */
      public double getX() {
         return x;
      }

      /**
       * @return the y
       */
      public double getY() {
         return y;
      }

      /**
       * @return the n
       */
      public double getN() {
         return n;
      }
   }

   /**
    * Creates new form Cyclotron
    */
   public Cyclotron(double radius, double current, double s, double object_mass, double magnetic_induction, double objects_per_second, double px_per_m) {
      initComponents();
      this.radius_m = radius;
      this.current = current;
      this.s = s;
      this.object_mass = object_mass;
      this.magnetic_induction = magnetic_induction;
      this.objects_per_second = objects_per_second;
      this.px_per_m = px_per_m;
      addFlyingObject();
   }

   public Cyclotron() {
      this(2, 100E-9, 1, 10, 0.000001, 1, 40);
   }

   /**
    * This method is called from within the constructor to initialize the form.
    * WARNING: Do NOT modify this code. The content of this method is always
    * regenerated by the Form Editor.
    */
   @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 400, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );
    }// </editor-fold>//GEN-END:initComponents
    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables

   @Override
   public void paintComponent(Graphics g) {
      g.clearRect(0, 0, getWidth(), getHeight());
      int width = (int) (px_per_m * (radius_m + s) * 2);
      int radius = (int) (px_per_m * radius_m);
      int s_px = (int) (px_per_m * s);
      g.setColor(Color.LIGHT_GRAY);
      g.fillArc((getWidth() - width) / 2 + (int) (s_px / 2.0), getHeight() / 2 - radius, radius * 2, radius * 2, 90, 180);
      g.fillArc((getWidth() - width) / 2 + (int) (3.0 / 2.0 * s_px), getHeight() / 2 - radius, radius * 2, radius * 2, 90, -180);
      g.setColor(Color.BLACK);
      g.drawString("Highest current pace = " + String.format("%4.3fm/s; radius = %4.3fm; n = %4.0f; E = %4.3fE-19J", max_v, max_r, max_n - 4, 0.5 * this.object_mass * Math.pow(max_v, 2) * 1.67), 5, 15);
      int middle_x = getWidth() / 2 - (int) (object_radius / 2);
      int middle_y = getHeight() / 2 - (int) (object_radius / 2);
      g.setColor(Color.red);
      try {
         for (FlyingObject flyingObject : flying_objects) {
            g.fillOval((int) ((flyingObject.getX() * px_per_m) + middle_x), (int) ((flyingObject.getY() * px_per_m) + middle_y), object_radius, object_radius);
         }
      } catch (Exception exp) {
      }
   }

   public void calc() {
      if (last_time == -1) {
         last_time = System.currentTimeMillis();
      }
      try {
         double add_time = (System.currentTimeMillis() - last_time) / 1000.0 * Math.pow(Math.abs(time_fac) + 1, Math.signum(time_fac));
         ArrayList<FlyingObject> del_list = new ArrayList<FlyingObject>();
         max_v = 0;
         max_n = 4;
         max_r = 0;
         for (FlyingObject flyingObject : flying_objects) {
            flyingObject.updateCoordinates(add_time);
            if (flyingObject.isOutOfZyklotron()) {
               del_list.add(flyingObject);
            }
         }
         flying_objects.removeAll(del_list);
      } catch (Exception exp) {
      }
      last_time = System.currentTimeMillis();
   }

   public void clear() {
      try {
         flying_objects.clear();
      } catch (Exception exp) {
      }
   }

   /**
    * @return the radius_m
    */
   public double getRadius_m() {
      return radius_m;
   }

   /**
    * @param radius_m the radius_m to set
    */
   public void setRadius_m(double radius_m) {
      this.radius_m = radius_m;
   }

   /**
    * @return the px_per_m
    */
   public double getPx_per_m() {
      return px_per_m;
   }

   /**
    * @param px_per_m the px_per_m to set
    */
   public void setPx_per_m(double px_per_m) {
      this.px_per_m = px_per_m;
   }

   /**
    * @return the s
    */
   public double getS() {
      return s;
   }

   /**
    * @param s the s to set
    */
   public void setS(double s) {
      this.s = s;
   }

   /**
    * @return the magnetic_induction
    */
   public double getMagnetic_induction() {
      return magnetic_induction;
   }

   /**
    * @param magnetic_induction the magnetic_induction to set
    */
   public void setMagnetic_induction(double magnetic_induction) {
      this.magnetic_induction = magnetic_induction;
   }

   /**
    * @return the current
    */
   public double getCurrent() {
      return current;
   }

   /**
    * @param current the current to set
    */
   public void setCurrent(double current) {
      this.current = current;
   }

   public void addFlyingObject() {
      flying_objects.add(new FlyingObject(object_mass));
   }

   public void setRadiusInMeter(double radius) {
      radius_m = radius;
   }
}
