/* //#region Assembly System.Numerics, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089
//// System.Numerics.dll
//#endregion
// using System;
using static System.Math;
using System.Collections;

 using  System.Numerics;

 using System.Collections.Generic;

namespace myVectors{
  public struct Vec2// : IEquatable<Vec2>, IFormattable
  {
    public double X;
    public double Y;

    public Vec2(double value) { X=value;Y=value;}
    public Vec2(double x, double y) {X=x;Y=y;}
//   public static Vec2 UnitX { get => new Vec2(1.0,0.0); }
//   public static Vec2 One { get => new Vec2(1.0,1.0); }
//   public static Vec2 UnitY { get => new Vec2(0.0,1.0); }
//   public static Vec2 Zero { get => new Vec2(0.0,0.0); }

  //  public static Vec2 Abs(Vec2 value);
  //  public static Vec2 Add(Vec2 left, Vec2 right) {return new Vec2(left.X+right.X,left.Y+right.Y);}
//       public static Vec2 Clamp(Vec2 value1, Vec2 min, Vec2 max);
//   public static double Distance(Vec2 value1, Vec2 value2);
//   public static double DistanceSquared(Vec2 value1, Vec2 value2);
//   public static Vec2 Divide(Vec2 left, Vec2 right);
//   public static Vec2 Divide(Vec2 left, double divisor);
//   public static double Dot(Vec2 value1, Vec2 value2);
//   public static Vec2 Lerp(Vec2 value1, Vec2 value2, double amount);
//   public static Vec2 Max(Vec2 value1, Vec2 value2);
//   public static Vec2 Min(Vec2 value1, Vec2 value2);
//   public static Vec2 Multiply(double left, Vec2 right);
//   public static Vec2 Multiply(Vec2 left, double right);
//   public static Vec2 Multiply(Vec2 left, Vec2 right);
//   public static Vec2 Negate(Vec2 value);
//   public static Vec2 Normalize(Vec2 value);
//   public static Vec2 Reflect(Vec2 vector, Vec2 normal);
//   public static Vec2 SquareRoot(Vec2 value);
//   public static Vec2 Subtract(Vec2 left, Vec2 right);
//   public static Vec2 Transform(Vec2 value, Quaternion rotation);
//   public static Vec2 Transform(Vec2 position, Matrix4x4 matrix);
//   public static Vec2 Transform(Vec2 position, Matrix3x2 matrix);
//   public static Vec2 TransformNormal(Vec2 normal, Matrix3x2 matrix);
//   public static Vec2 TransformNormal(Vec2 normal, Matrix4x4 matrix);
//   public void CopyTo(double[] array, int index);
//   public void CopyTo(double[] array);
//   public bool Equals(Vec2 other);
//   public override bool Equals(object obj);
//   public override int GetHashCode();
//   public double Length();
//   public double LengthSquared();
//   public string ToString(string format, IFormatProvider formatProvider);
//   public override string ToString();
//   public string ToString(string format);
//
  public static Vec2 operator +(Vec2 left, Vec2 right){return new Vec2(left.X+right.X,left.Y+right.Y);}
  public static Vec2 operator -(Vec2 value){return new Vec2(-value.X,-value.Y);}
  public static Vec2 operator -(Vec2 left, Vec2 right){return new Vec2(left.X-right.X,left.Y-right.Y);}
  public static Vec2 operator *(Vec2 left, double right){return new Vec2(left.X*right,left.Y*right);}
  public static Vec2 operator *(Vec2 left, Vec2 right){return new Vec2(left.X*right.X,left.Y*right.Y);}
  public static Vec2 operator *(double left, Vec2 right){return new Vec2(right.X*left,right.Y*left);}
  public static Vec2 operator /(Vec2 left, Vec2 right){return new Vec2(left.X/right.X,left.Y/right.Y);}
  public static Vec2 operator /(Vec2 value1, double value2){return new Vec2(value1.X/value2,value1.Y/value2);}
// public static bool operator ==(Vec2 left, Vec2 right);
// public static bool operator !=(Vec2 left, Vec2 right);
  }
  public struct Vec3// : IEquatable<Vec2>, IFormattable
  {
    public double X;
    public double Y;
    public double Z;
    public Vec3(double value) { X=value;Y=value;Z=value;}
    public Vec3(double x, double y, double z) {X=x;Y=y;Z=z;}

public static Vec3 operator +(Vec3 left, Vec3 right){return new Vec3(left.X+right.X,left.Y+right.Y,left.Z+right.Z);}
public static Vec3 operator -(Vec3 value){return new Vec3(-value.X,-value.Y,-value.Z);}
public static Vec3 operator -(Vec3 left, Vec3 right){return new Vec3(left.X-right.X,left.Y-right.Y,left.Z-right.Z);}
public static Vec3 operator *(Vec3 left, double right){return new Vec3(left.X*right,left.Y*right,left.Z*right);}
public static Vec3 operator *(Vec3 left, Vec3 right){return new Vec3(left.X*right.X,left.Y*right.Y,left.Z*right.Z);}
public static Vec3 operator *(double left, Vec3 right){return new Vec3(right.X*left,right.Y*left,right.Z*left);}
public static Vec3 operator /(Vec3 left, Vec3 right){return new Vec3(left.X/right.X,left.Y/right.Y,left.Z/right.Z);}
public static Vec3 operator /(Vec3 value1, double value2){return new Vec3(value1.X/value2,value1.Y/value2,value1.Z/value2);}
// public static bool operator ==(Vec2 left, Vec2 right);
// public static bool operator !=(Vec2 left, Vec2 right);
  }
} */