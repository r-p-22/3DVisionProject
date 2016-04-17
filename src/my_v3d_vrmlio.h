#ifndef V3D_VRMLIO_H
#define V3D_VRMLIO_H

#include "camera.h"
#include <iostream>
//#include <fstream>
#include <string>
#include <cstdio>

    template<typename U>
    void writeAppearanceToVRMLStream(ofstream & os,U const& color255, float trans=0.5);

   template <typename T,typename U>
   void writeColoredPointsToVRML(std::vector<T> const& points, std::vector<U> const& colors255,
                                 char const * filename, bool append = false)
   {
      using namespace std;
      std::ofstream os;

      os.open(filename, append ? ios::app : ios::out);
      os.precision(15);

      if (!append) os <<"#VRML V2.0 utf8" << endl;
      os << "Shape {" << endl;

      os << " appearance Appearance {" << endl;
      os << "   material Material {" << endl;
      os << "        ambientIntensity 0.1" << endl;
      os << "        diffuseColor 1.0 0.7 0.2" << endl;
      os << " } } " << endl;
      os << "  geometry PointSet {" << endl;

      // write point coordinates
      os << "  coord Coordinate {" << endl;
      os << "  point [" << endl;
      for (size_t i = 0; i < points.size(); ++i)
         os << points[i][0] <<" " << points[i][1] << " " << points[i][2] << "," << endl;
      os << " ]\n } # end points"<<endl;

      os << " color Color { color [" << endl;
      for (size_t i = 0; i < points.size(); ++i)
      {
         double r = ((double) colors255[i][0]) / 255.0;
         double g = ((double) colors255[i][1]) / 255.0;
         double b = ((double) colors255[i][2]) / 255.0;
         os << min(max(r, 0.0), 1.0) << " " << min(max(g, 0.0), 1.0) << " " << min(max(b, 0.0), 1.0) << " ," << endl;
      }
      os << " ] } # end color" << endl;
      os << " }" << endl;
      os << "}" << endl;
   } // end writeColoredPointsToVRML()

   template <typename U>
   void writeCameraFrustumToVRML(CameraMatrix const& cam, int width, int height, double scale,
                                 U const& color255, const char *filename, bool append)
   {
      using namespace std;

      ofstream os;

      os.open(filename, append ? ios::app : ios::out);
      os.precision(15);

      double r = ((double) color255[0]) / 255.0;
      double g = ((double) color255[1]) / 255.0;
      double b = ((double) color255[2]) / 255.0;

      if (!append) os <<"#VRML V2.0 utf8" << endl;
      os << " Shape {" << endl;
      os << " appearance Appearance {" << endl;
      os << "  material Material { emissiveColor 0.2 0.2 0.2" << endl;
      os << "  diffuseColor " << r << " " << g << " " << b << endl;
      os << "  }" << endl;
      os << " }" << endl;

      os <<" geometry IndexedFaceSet {" << endl;
      os <<" solid FALSE " << endl;
      os <<" ccw TRUE " << endl;
      os <<" coord Coordinate {" << endl;
      os <<" point [" << endl;

      Vector3d const center = cam.cameraCenter();
      Vector3d const ll = cam.unprojectPixel(Vector2f(0.0, double(height)), scale);
      Vector3d const lr = cam.unprojectPixel(Vector2f(double(width), double(height)), scale);
      Vector3d const ul = cam.unprojectPixel(Vector2f(0.0, 0.0), scale);
      Vector3d const ur = cam.unprojectPixel(Vector2f(double(width), 0.0), scale);

      // write point coordinates
      os << ll[0] << " " << ll[1] << " " << ll[2] << "," << endl;
      os << lr[0] << " " << lr[1] << " " << lr[2] << "," << endl;
      os << ul[0] << " " << ul[1] << " " << ul[2] << "," << endl;
      os << ur[0] << " " << ur[1] << " " << ur[2] << "," << endl;
      os << center[0] << " " << center[1] << " " << center[2] << "," << endl;

      os << " ]\n } "<< endl;

      os << "  coordIndex [\n";
      os << "  0,1,2,3,-1" << endl;
      os << "  0,1,4,-1" << endl;
      os << "  1,3,4,-1" << endl;
      os << "  3,2,4,-1" << endl;
      os << "  2,0,4,-1\n  ]" << endl;
      os << " }\n}" << endl;
   } // end writeCameraFrustumToVRML()

   template <typename T>
   void writeMatchesToVRML(std::vector<T> const& points, char const * filename, bool append = false)
   {
       using namespace std;
       std::ofstream os;

       os.open(filename, append ? ios::app : ios::out);
       os.precision(15);

       os << "#VRML V2.0 utf8" << endl;
       os << " Shape {" << endl;

       os << " appearance Appearance {" << endl;
       os << " material Material {" << endl;
       os << "     ambientIntensity 0.1" << endl;
       os << "     diffuseColor 1.0 0.0 0.0" << endl;
       os << "     transparency      0.5" << endl;
       os << " } } " << endl;
       os << "  geometry IndexedLineSet {" << endl;

       // write lines between matches coordinates
       os <<"  coord Coordinate {" << endl;
       os <<"  point [" << endl;
       for (size_t i = 0; i < points.size(); ++i)
           os << points[i][0] <<" " << points[i][1] << " " << points[i][2] << "," << endl;
       os << " ]\n } # end coord"<<endl;

       os <<"  coordIndex ["<<endl;
       for(size_t i=0; i<points.size()-1; i+=2)
           os <<i<<", "<<i+1<<", -1,"<<endl;
       os <<"  ] #end coordIndex"<<endl;

       os << " } # end geometry" << endl;
       os << "} # end shape" << endl;
   }

   template <typename T,typename U>
   void writePointsToVRML(std::vector<T> const& points, U const & color255,char const * filename, bool append = false)
   {
	   using namespace std;
	         std::ofstream os;

	         os.open(filename, append ? ios::app : ios::out);
	         os.precision(15);

	         if (!append) os <<"#VRML V2.0 utf8" << endl;
	         os << "Shape {" << endl;

	         os << " appearance Appearance {" << endl;
	         os << "   material Material {" << endl;
	         os << "        ambientIntensity 0.1" << endl;
	         os << "        diffuseColor 1.0 0.7 0.2" << endl;
	         os << " } } " << endl;
	         os << "  geometry PointSet {" << endl;

	         // write point coordinates
	         os << "  coord Coordinate {" << endl;
	         os << "  point [" << endl;
	         for (size_t i = 0; i < points.size(); ++i)
	            os << points[i][0] <<" " << points[i][1] << " " << points[i][2] << "," << endl;
	         os << " ]\n } # end points"<<endl;

	         os << " color Color { color [" << endl;
	         for (size_t i = 0; i < points.size(); ++i)
	         {
	            double r = ((double) color255[0]) / 255.0;
	            double g = ((double) color255[1]) / 255.0;
	            double b = ((double) color255[2]) / 255.0;
	            os << min(max(r, 0.0), 1.0) << " " << min(max(g, 0.0), 1.0) << " " << min(max(b, 0.0), 1.0) << " ," << endl;
	         }
	         os << " ] } # end color" << endl;
	         os << " }" << endl;
	         os << "}" << endl;
   }

   template <typename T>
   void writePointsToVRML(std::vector<T> const& points, char const * filename, bool append = false)
   {
      Vector3f color255(0.,0.,255.0);
      writePointsToVRML(points, color255, filename, append);
   }

   template <typename U>
   void writeCameraFrustumToVRMLStream(ofstream & os,CameraMatrix const& cam,  double scale,
                                  U const& color255)
    {
       using namespace std;

       os.precision(3);
       os << " Shape {" << endl;
       writeAppearanceToVRMLStream(os,color255);

       os.precision(15);
       os <<" geometry IndexedFaceSet {" << endl;
       os <<" solid FALSE " << endl;
       os <<" ccw TRUE " << endl;
       os <<" coord Coordinate {" << endl;
       os <<" point [" << endl;
       Vector3d const center = cam.cameraCenter();
       Vector3d ll = cam.unprojectPixel(Vector2f(0.0, (double)cam.getHeight()), scale);
       Vector3d lr = cam.unprojectPixel(Vector2f((double)cam.getWidth(), (double)cam.getHeight()), scale);
       Vector3d ul = cam.unprojectPixel(Vector2f(0.0, 0.0), scale);
       Vector3d ur = cam.unprojectPixel(Vector2f((double)cam.getWidth(), 0.0), scale);
       Matrix3d I=Matrix3d::Identity();
       if(cam.getIntrinsic()==I || cam.getHeight()==0){//If normalized Camera, assume 4/3 AspectRatio for Frustrum
          ll = cam.unprojectNormalizedPixel(Vector2f(-0.5,  0.375), scale);
          lr = cam.unprojectNormalizedPixel(Vector2f( 0.5,  0.375), scale);
          ul = cam.unprojectNormalizedPixel(Vector2f(-0.5, -0.375), scale);
          ur = cam.unprojectNormalizedPixel(Vector2f( 0.5, -0.375), scale);
          os<<"#normalized camera assumed."<<endl;
       }
       // write point coordinates
       os << ll[0] << " " << ll[1] << " " << ll[2] << "," << endl;
       os << lr[0] << " " << lr[1] << " " << lr[2] << "," << endl;
       os << ul[0] << " " << ul[1] << " " << ul[2] << "," << endl;
       os << ur[0] << " " << ur[1] << " " << ur[2] << "," << endl;
       os << center[0] << " " << center[1] << " " << center[2] << "," << endl;

       os << " ]\n } "<< endl;

       os << "  coordIndex [\n";
       //os << "  0,1,2,3,-1" << endl;
       os << "  0,1,4,-1" << endl;
       os << "  1,3,4,-1" << endl;
       os << "  3,2,4,-1" << endl;
       os << "  2,0,4,-1\n  ]" << endl;
       os << " }\n}" << endl;
    } // end writeCameraFrustumToVRML()

   template <typename U>
   void writeCameraImageToVRMLStream(ofstream & os,CameraMatrix const& cam,  double scale,
                                  U const & imageName)
    {
       using namespace std;

       os.precision(3);
       os << "Shape {" << endl;
       os << " appearance Appearance {"<< endl;
       os << "   texture ImageTexture { url \""<<imageName<<"\"}\n }"<< endl;

       os.precision(15);
       os <<" geometry IndexedFaceSet {" << endl;
       //os <<"  texCoordIndex [0 1 3 -1 3 1 2 -1 ]" << endl; // We might need to flip thos if image is upside down
       os <<"  texCoordIndex [0 3 1 -1 1 3 2]" << endl; // We might need to flip thos if image is upside down
       os <<"  texCoord TextureCoordinate {point [0 0 1 0 1 1 0 1 ]}" << endl;

       Vector3d const center = cam.cameraCenter();
       Vector3d ll = cam.unprojectPixel(Vector2f(0.0, (double)cam.getHeight()), scale);
       Vector3d lr = cam.unprojectPixel(Vector2f((double)cam.getWidth(), (double)cam.getHeight()), scale);
       Vector3d ul = cam.unprojectPixel(Vector2f(0.0, 0.0), scale);
       Vector3d ur = cam.unprojectPixel(Vector2f((double)cam.getWidth(), 0.0), scale);
       Matrix3d I=Matrix3d::Identity();
       if(cam.getIntrinsic()==I || cam.getHeight()==0){//If normalized Camera, assume 4/3 AspectRatio for Frustrum
          ll = cam.unprojectNormalizedPixel(Vector2f(-0.5,  0.375), scale);
          lr = cam.unprojectNormalizedPixel(Vector2f( 0.5,  0.375), scale);
          ul = cam.unprojectNormalizedPixel(Vector2f(-0.5, -0.375), scale);
          ur = cam.unprojectNormalizedPixel(Vector2f( 0.5, -0.375), scale);
          os<<"#normalized camera assumed."<<endl;
       }
       // write point coordinates
       os <<"  coord Coordinate {" << endl;
       os <<"   point [" << endl;
       os << ll[0] << " " << ll[1] << " " << ll[2] << "," << endl;
       os << lr[0] << " " << lr[1] << " " << lr[2] << "," << endl;
       os << ul[0] << " " << ul[1] << " " << ul[2] << "," << endl;
       os << ur[0] << " " << ur[1] << " " << ur[2] << "," << endl;
       os << " ]\n  } "<< endl;

       os << "  coordIndex [0 2 1 -1 1 2 3 -1]" << endl;
       os << " }\n}" << endl;
    } // end writeCameraImageToVRML()

   template <typename U>
    void writeCameraFrustumToVRML(CameraMatrix const& cam, double scale,
                                   U const& color255, const char *filename, bool append)
     {
        using namespace std;
        ofstream os;
        os.open(filename, append ? ios::app : ios::out);
        os.precision(3);
        if (!append) os <<"#VRML V2.0 utf8" << endl;
        writeCameraFrustumToVRMLStream(os,cam,  scale, color255);
     } // end writeCameraFrustumToVRML()

    template <typename U>
     void writeCameraFrustumAndImageToVRML(CameraMatrix const& cam, double scale,
                                    U const& color255, std::string & imageName,const char *filename, bool append)
      {
         using namespace std;
         ofstream os;
         os.open(filename, append ? ios::app : ios::out);
         os.precision(3);
         if (!append) os <<"#VRML V2.0 utf8" << endl;
         writeCameraFrustumToVRMLStream(os,cam,  scale, color255);
         writeCameraImageToVRMLStream(os,cam,scale*2,imageName);
      } // end writeCameraFrustumToVRML()

    template<typename U>
    void writeAppearanceToVRMLStream(ofstream & os,U const& color255, float transp=0.5)
    {
       double r = ((double) color255[0]) / 255.0;
       double g = ((double) color255[1]) / 255.0;
       double b = ((double) color255[2]) / 255.0;
       os << " appearance Appearance {" << endl;
       os << "  material Material { ambientIntensity 0.1 " << endl;
       os << "  diffuseColor " << r << " " << g << " " << b << endl;
       //os << "  transparency  "<< transp << endl;
       os << "  }" << endl;
       os << " }" << endl;
    }

    template <typename U,typename T,typename V>
    void write2DnormalizedVs3DMatchesToVRMLStream(ofstream & os,CameraMatrix const& cam,  std::vector< Vector2f > const& points2D,std::vector< Vector3d > const& points3D,double scale,
                                   vector< Vector3f > const& colors255)
     {
        using namespace std;
        os.precision(15);
        assert(points3D.size()==points2D.size());
        assert(colors255.size()==points2D.size());
        os << " Shape {" << endl;
        Vector3i color255(255,255,255);
        writeAppearanceToVRMLStream(os,color255);
        os << "  geometry IndexedFaceSet {" << endl;

        // write lines between matches coordinates
        os <<"  coord Coordinate {" << endl;
        os <<"  point [" << endl;
        os << cam.cameraCenter()[0] <<" " << cam.cameraCenter()[1] << " " << cam.cameraCenter()[2] << "," << endl;
        for (size_t i = 0; i < points2D.size(); ++i){
           Vector3d p2=cam.unprojectPixel(points2D[i],scale);
           os << p2[0] <<" " << p2[1] << " " << p2[2] << "," << endl;
           os << p2[0] <<" " << p2[1] << " " << p2[2] << "," << endl;
           os << points3D[i][0] <<" " << points3D[i][1] << " " << points3D[i][2] << "," << endl;
        }
        os << " ]\n } # end coord"<<endl;

        os <<"  coordIndex ["<<endl;
        for(size_t i=0; i<points2D.size(); i++){
           os <<0<<", "<<3*i+1<<", "<<3*i+2<<", -1,"<<endl;
           os <<3*i+1<<", "<<3*i+2<<", "<<3*i+3<<", -1,"<<endl;
        }
        os <<"  ] #end coordIndex"<<endl;

        os << " colorPerVertex TRUE" << endl;
        os << " color Color { color [" << endl;
        os << 0 << " " << 0 << " " << 0 << " ," << endl;
        for (size_t i = 0; i < points2D.size(); ++i)
        {
           double r = min(max(((double) colors255[i][0]) / 255.0, 0.0), 1.0);
           double g = min(max(((double) colors255[i][1]) / 255.0, 0.0), 1.0);
           double b = min(max(((double) colors255[i][2]) / 255.0, 0.0), 1.0);
           os << r << " " << g << " " << b << " ," << endl;
           os << r << " " << g << " " << b << " ," << endl;
           os << r << " " << g << " " << b << " ," << endl;
        }
        os << " ] } # end Color" << endl;
        os << " } # end geometry" << endl;
        os << "} # end shape" << endl;

     } // end writeCamera3D2DMatchesToVRMLStream()
    template <typename U,typename T,typename V>
    void write2DnormalizedVs3DMatchesToVRMLStream(ofstream & os,CameraMatrix const& cam,  std::vector<Vector2f > const& points2D,std::vector<Vector3d > const& points3D,double scale,
                                   Vector3f const& color255)
     {
        using namespace std;
        os.precision(15);
        CameraMatrix normalizedCam(cam);
        Matrix3d K=Matrix3d::Identity();normalizedCam.setIntrinsic(K);

        os << " Shape {" << endl;
        writeAppearanceToVRMLStream(os,color255);
        os << "  geometry IndexedFaceSet {" << endl;

        // write lines between matches coordinates
        os <<"  coord Coordinate {" << endl;
        os <<"  point [" << endl;
        os << normalizedCam.cameraCenter()[0] <<" " << normalizedCam.cameraCenter()[1] << " " << normalizedCam.cameraCenter()[2] << "," << endl;
        for (size_t i = 0; i < points2D.size(); ++i){
           Vector3d p2=normalizedCam.unprojectPixel(points2D[i],scale);
           os << p2[0] <<" " << p2[1] << " " << p2[2] << "," << endl;
           os << p2[0] <<" " << p2[1] << " " << p2[2] << "," << endl;
           os << points3D[i][0] <<" " << points3D[i][1] << " " << points3D[i][2] << "," << endl;
        }
        os << " ]\n } # end coord"<<endl;

        os <<"  coordIndex ["<<endl;
        for(size_t i=0; i<points2D.size(); i++){
           os <<0<<", "<<3*i+1<<", "<<3*i+2<<", -1,"<<endl;
           os <<3*i+1<<", "<<3*i+2<<", "<<3*i+3<<", -1,"<<endl;
        }
        os <<"  ] #end coordIndex"<<endl;

        os << " } # end geometry" << endl;
        os << "} # end shape" << endl;

     } // end writeCamera3D2DMatchesToVRMLStream()

   template <typename U,typename T,typename V>
   void write2DnormalizedVs3DMatchesToVRML(CameraMatrix const& cam,  std::vector<Vector2f > const& points2D,std::vector<Vector3d > const& points3D,double scale,
                                  U const& color255,const char *filename, bool append)
   {
      using namespace std;
      ofstream os;
      os.open(filename, append ? ios::app : ios::out);
      os.precision(3);
      if (!append) os <<"#VRML V2.0 utf8" << endl;
      write2DnormalizedVs3DMatchesToVRMLStream(os, cam, points2D, points3D, scale, color255);
   }

   template<typename T>

   inline void
   writeQuantilePointsToVRML(std::vector<T> const &goodXs, char const * wrlName, float quantile, bool append=false)
   {
      using namespace std;

      vector<Vector3d> goodXs2;

      vector<float> norms(goodXs.size());

      for (size_t i = 0; i < goodXs.size(); ++i) norms[i] = goodXs[i].norm();
      std::sort(norms.begin(), norms.end());
      float distThr = norms[int(norms.size() * quantile)];
      //cout << "90% quantile distance: " << distThr << endl;

      for (size_t i = 0; i < goodXs.size(); ++i)
      {
         Vector3d const& X = goodXs[i];
         if (X.norm() < 3*distThr) goodXs2.push_back(X);
      }
      writePointsToVRML(goodXs2,  wrlName, append);
   } // end writeQuantilePointsToVRML()

   template <typename T,typename U>

   inline void
   writeQuantilePointsToVRML(std::vector<T> const &goodXs, U &color, char const * wrlName, float quantile, bool append=false)
   {
      using namespace std;

      vector<Vector3d> goodXs2;

      vector<float> norms(goodXs.size());

      for (size_t i = 0; i < goodXs.size(); ++i) norms[i] = goodXs[i].squaredNorm();
      std::sort(norms.begin(), norms.end());
      float distThr = norms[int(norms.size() * quantile)];
      //cout << "90% quantile distance: " << distThr << endl;

      for (size_t i = 0; i < goodXs.size(); ++i)
      {
         Vector3d const& X = goodXs[i];
         if (X.norm() < 3*distThr) goodXs2.push_back(X);
      }

      writePointsToVRML(goodXs2,color,  wrlName, append);
   } // end writeQuantilePointsToVRML()



   //>-------------------------------------------------------------------------------->>
   //>-------------------------------------------------------------------------------->>
   //>-------------------------------------------------------------------------------->>
   //>-------------------------------------------------------------------------------->>
   //>-------------------------------------------------------------------------------->>
   //>-------------------------------------------------------------------------------->>

   template<typename T>
   void writeGroupsToVRML(std::vector< std::vector<T> > groups,char const* wrlName,
		   float quantile)
      {

	   int numgroups = groups.size();

	   float q = 255./numgroups;

	   for (int i=0; i < numgroups; i++){

		   Eigen::Vector3f color255(q*i, 255. - q*i, 0.);
		   writeQuantilePointsToVRML(groups[i], color255,wrlName, quantile, i>0);

	   }
      }


   template<typename T>
   void writePlanesToVRML(std::vector<T> const points3d, std::vector<Eigen::Vector4d> planes,
  		char const* wrlName, float quantile, bool append=false)
   {
     using namespace std;

     vector<Eigen::Vector3d> goodXs2;

	 vector<float> norms(points3d.size());

	 for (size_t i = 0; i < points3d.size(); ++i) norms[i] = points3d[i].squaredNorm();
	 std::sort(norms.begin(), norms.end());
	 float distThr = norms[int(norms.size() * quantile)];
	 //cout << "90% quantile distance: " << distThr << endl;

	 for (size_t i = 0; i < points3d.size(); ++i)
	 {
		 Eigen::Vector3d const& X = points3d[i];
		if (X.norm() < 3*distThr) goodXs2.push_back(X);
	 }

	 //forget about points3d, goodXs2 is the point set

	 //open file stream
	 Vector3f color255(255.0,0.,255.0);
	 std::ofstream os;
     os.open(wrlName, append ? ios::app : ios::out);

	 //Get the point cloud bounding points
	 //operate on Homogeneous coordinates
	 Eigen::Vector4d bboxMin; bboxMin.head(3) = (goodXs2[0]); bboxMin(3) = 1.;
	 Eigen::Vector4d bboxMax; bboxMax.head(3) = (goodXs2[0]); bboxMax(3) = 1.;
	 for (int j = 1; j < goodXs2.size(); ++j){
		 bboxMin[0] = std::min(points3d[j][0], bboxMin[0]);
		 bboxMin[1] = std::min(points3d[j][1], bboxMin[1]);
		 bboxMin[2] = std::min(points3d[j][2], bboxMin[2]);
		 bboxMax[0] = std::max(points3d[j][0], bboxMax[0]);
		 bboxMax[1] = std::max(points3d[j][1], bboxMax[1]);
		 bboxMax[2] = std::max(points3d[j][2], bboxMax[2]);
	  }

	 cout << "planes SIZE: " << planes.size() << endl;

	 Eigen::Vector4d testplane;
	 for (int pli=0;pli<planes.size();pli++){
		 testplane = planes.at(pli);
		 std::vector<Eigen::Vector3d> planepoints;

		 Eigen::Vector3d v;
		 v << bboxMin[0], bboxMin[1], 0;
		 v(2) = (-testplane(3)-v(0)*testplane(0)-v(1)*testplane(1))/testplane(2);
		 planepoints.push_back(v);

		 //cout << "min==" << endl;
		 //cout << v << endl;

		 v << bboxMax[0], bboxMin[1], 0;
		 v(2) = (-testplane(3)-v(0)*testplane(0)-v(1)*testplane(1))/testplane(2);
		 planepoints.push_back(v);

		 v << bboxMin[0], bboxMax[1], 0;
		 v(2) = (-testplane(3)-v(0)*testplane(0)-v(1)*testplane(1))/testplane(2);
		 planepoints.push_back(v);

		 v << bboxMax[0], bboxMax[1], 0;
		 v(2) = (-testplane(3)-v(0)*testplane(0)-v(1)*testplane(1))/testplane(2);
		 planepoints.push_back(v);

		 //cout << "max==" << endl;
		 //cout << v << endl;
		 //write the plane to VRML file


		 os.precision(3);
		 os << " Shape {" << endl;
		 writeAppearanceToVRMLStream(os,color255);

		os.precision(15);
		os <<" geometry IndexedFaceSet {" << endl;
		os <<" solid FALSE " << endl;
		os <<" ccw TRUE " << endl;
		os <<" coord Coordinate {" << endl;
		os <<" point [" << endl;

		// write point coordinates
		Eigen::Vector3d ll = planepoints[0];
		Eigen::Vector3d lr = planepoints[1];
		Eigen::Vector3d ul = planepoints[2];
		Eigen::Vector3d ur = planepoints[3];
		os << ll[0] << " " << ll[1] << " " << ll[2] << "," << endl;
		os << lr[0] << " " << lr[1] << " " << lr[2] << "," << endl;
		os << ul[0] << " " << ul[1] << " " << ul[2] << "," << endl;
		os << ur[0] << " " << ur[1] << " " << ur[2] << "," << endl;

		os << " ]\n } # end coord"<<endl;

		os << "  coordIndex [\n";
		//os << "  0,1,2,3,-1" << endl;
		os << "  0,1,2,-1" << endl;
		os << "  1,2,3,-1" << endl;
		//os << "  2,3,0,-1" << endl;
		os << "  0,2,3,-1\n  ]" << endl;
		os << " }\n}" << endl;

	 }
	 os.close();
   }

   using namespace std;

   void writeLatticeToVRML(Vector4d plane,vector<Eigen::Vector3d> basisVectors, vector<Eigen::Vector3d> Bounds,
		   char const * filename, bool append=true){

		ofstream os;
		os.open(filename, append ? ios::app : ios::out);
		os.precision(4);
		if (!append) os <<"#VRML V2.0 utf8" << endl;

		//Assume bound points ordering:
		//1st lowerLeft: LL
		//2nd topRight: TR
		Vector3d LL = Bounds[0];
		Vector3d TR = Bounds[1];
		Vector3d basis1 = Bounds[0];
		Vector3d basis2 = Bounds[1];

		//Define topLeft and lowerRight points
		double cos_phi1 = (TR-LL).dot(basis1)/( sqrt((TR-LL).squaredNorm())*sqrt(basis1.squaredNorm()));
		Vector3d B1; B1 = sqrt((TR-LL).squaredNorm())*cos_phi1*basis1/sqrt(basis1.squaredNorm()) + LL;

		double cos_phi2 = (TR-LL).dot(basis2)/( sqrt((TR-LL).squaredNorm())*sqrt(basis2.squaredNorm()));
		Vector3d B2; B2 = sqrt((TR-LL).squaredNorm())*cos_phi2*basis2/sqrt(basis2.squaredNorm()) + LL;

		//number of lattices in the 2 axes
		int k1 = round( sqrt((B1-LL).squaredNorm())/sqrt(basis1.squaredNorm()) );
		int k2 = round( sqrt((B2-LL).squaredNorm())/sqrt(basis2.squaredNorm()) );

		cout << k1 << endl;
		cout << k2 << endl;
		os << " Shape {" << endl;

		os << " appearance Appearance {" << endl;
		os << " material Material {" << endl;
		os << "     ambientIntensity 0.1" << endl;
		os << "     diffuseColor 1.0 0.0 0.0" << endl;
		os << "     transparency      0.5" << endl;
		os << " } } " << endl;
		os << "  geometry IndexedLineSet {" << endl;

		// write lines between matches coordinates
		os <<"  coord Coordinate {" << endl;
		os <<"  point [" << endl;

		for (int k=0; k<=k2; k++){
			Vector3d pa; pa = LL + k*basis2;
			Vector3d pb; pb = B1 + k*basis2;

			os << pa[0] <<" " << pa[1] << " " << pa[2] << "," << endl;
			os << pb[0] <<" " << pb[1] << " " << pb[2] << "," << endl;

		}
		for (int k=0; k<=k1; k++){
			Vector3d pa; pa = LL + k*basis1;
			Vector3d pb; pb = B2 + k*basis1;

			os << pa[0] <<" " << pa[1] << " " << pa[2] << "," << endl;
			os << pb[0] <<" " << pb[1] << " " << pb[2] << "," << endl;
		}

		os << " ]\n } # end coord"<<endl;

		os << " color Color { color [" << endl;
		for(int i=0; i<1; i++)//i<2*(k1+k2+2); i+=2)
		{
			double r = 1.0;
			double g = 0.0;
			double b = 1.0;
			os << min(max(r, 0.0), 1.0) << " " << min(max(g, 0.0), 1.0) << " " << min(max(b, 0.0), 1.0) << " ," << endl;
		}
		os << " ] } # end color" << endl;


		os <<"  coordIndex ["<<endl;
		for(int i=0; i<2*(k1+k2+2); i+=2)
		   os <<i<<", "<<i+1<<", -1,"<<endl;
		os <<"  ] #end coordIndex"<<endl;

		os <<"  colorIndex ["<<endl;
			os <<"0 ";
			for(int i=2; i<2*(k1+k2+2); i+=2)
				os <<",0 ";
		os << endl <<"  ] #end colorIndex"<<endl;

		os << " } # end geometry" << endl;
		os << "} # end shape" << endl;

}



#endif // V3D_VRMLIO_H
