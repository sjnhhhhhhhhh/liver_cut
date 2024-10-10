#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>
#include <vtkPointPicker.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkProperty.h>
#include <vtkSTLReader.h>
#include <vtkCylinderSource.h>
#include <vtkCubeAxesActor.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAutoInit.h>
#include <vtkPlaneSource.h>
#include <vtkDataSetMapper.h>
#include <vtkSTLWriter.h>
#include <vtkGeometryFilter.h>
#include <vtkTriangle.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkContourFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkClipPolyData.h>
#include <vtkSampleFunction.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkLineSource.h>
#include <vtkConvexHull2D.h>
#include <vtkDelaunay3D.h>
#include <vtkGeometryFilter.h>
#include <vtkConvexPointSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCommand.h>
#include <vtkCamera.h>
#include <chrono>
#include <Eigen/Dense>
#include <utility>  // For std::move

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include "liver_cut.h"
#include "cut_tail.h"
#include "cut_right_lobe.h"
#include "cut_left_lobe.h"
#include "cut_rb3b.h"
#include "cut_rb2b.h"
#include "cut_rf_up.h"
#include "cut_rb1b_rf_down.h"
#include "cut_left_inside.h"
#include "cut_left_ou_od.h"

using namespace std;


int main() {
    // 读取STL文件
    vtkSmartPointer<vtkSTLReader> reader1 = vtkSmartPointer<vtkSTLReader>::New();
    vtkSmartPointer<vtkSTLReader> reader2 = vtkSmartPointer<vtkSTLReader>::New();
    vtkSmartPointer<vtkSTLReader> reader3 = vtkSmartPointer<vtkSTLReader>::New();
    vtkSmartPointer<vtkSTLReader> reader4 = vtkSmartPointer<vtkSTLReader>::New();
    vtkSmartPointer<vtkSTLReader> reader5 = vtkSmartPointer<vtkSTLReader>::New();

    reader1->SetFileName("C:/code/liver_cut/data/rb1v.stl");
    reader2->SetFileName("C:/code/liver_cut/data/rb2v.stl");
    reader3->SetFileName("C:/code/liver_cut/data/rb3v.stl");
    reader4->SetFileName("C:/code/liver_cut/data/rf_up_v.stl");
    reader5->SetFileName("C:/code/liver_cut/data/rf_down_v.stl");
    reader1->Update();
    reader2->Update();
    reader3->Update();
    reader4->Update();
    reader5->Update();

    vtkSmartPointer<vtkSTLReader> reader6 = vtkSmartPointer<vtkSTLReader>::New();
    reader6->SetFileName("C:/code/liver_cut/data/left_inside_v.stl");
    reader6->Update();
    auto left_inside_v = reader6->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader7 = vtkSmartPointer<vtkSTLReader>::New();
    reader7->SetFileName("C:/code/liver_cut/data/left_ou_v.stl");
    reader7->Update();
    auto left_ou_v = reader7->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader8 = vtkSmartPointer<vtkSTLReader>::New();
    reader8->SetFileName("C:/code/liver_cut/data/left_od_v.stl");
    reader8->Update();
    auto left_od_v = reader8->GetOutput();

    vtkSmartPointer<vtkSTLReader> reader9 = vtkSmartPointer<vtkSTLReader>::New();
    reader9->SetFileName("C:/code/liver_cut/data/tail_v.stl");
    reader9->Update();
    auto tail = reader9->GetOutput();
    



   /*vtkSmartPointer<vtkSTLReader> reader10 = vtkSmartPointer<vtkSTLReader>::New();
    reader10->SetFileName("C:/code/liver_cut/data/liver.stl");
    reader10->Update();
    auto liver = reader10->GetOutput();
    double mbounds_liver[6];
    liver->GetBounds(mbounds_liver);*/
    vtkSmartPointer<vtkPolyData> liver = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkSTLReader> reader10 = vtkSmartPointer<vtkSTLReader>::New();
        reader10->SetFileName("C:/code/liver_cut/data/1.stl");
        reader10->Update();
        liver = mergePolyData2(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/2.stl");
        reader10->Update();
        liver = mergePolyData2(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/3.stl");
        reader10->Update();
        liver = mergePolyData2(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/4a.stl");
        reader10->Update();
        liver = mergePolyData2(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/4b.stl");
        reader10->Update();
        liver = mergePolyData2(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/5.stl");
        reader10->Update();
        liver = mergePolyData2(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/6.stl");
        reader10->Update();
        liver = mergePolyData2(liver,reader10->GetOutput());
        
        reader10->SetFileName("C:/code/liver_cut/data/7.stl");
        reader10->Update();
        liver = mergePolyData2(liver,reader10->GetOutput());

        reader10->SetFileName("C:/code/liver_cut/data/8.stl");
        reader10->Update();
        liver = mergePolyData2(liver,reader10->GetOutput());   

        // 使用 vtkCleanPolyData 来移除重复的点和单元，实测做不到
        auto cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
        cleanFilter->SetInputData(liver);
        cleanFilter->Update();

        auto liver_cl = vtkSmartPointer<vtkPolyData>::New();
        liver_cl->ShallowCopy(cleanFilter->GetOutput());

        double mbounds_liver[6];
        liver_cl->GetBounds(mbounds_liver);

        std::cout<<liver_cl->GetNumberOfPoints()<<"\n";
        vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
        stlWriter->SetFileName("C:/code/liver_cut/new_stl/liver_clean.stl");
        stlWriter->SetInputData(liver_cl);
        stlWriter->Write();

    // 提取读取的数据
    auto rb1v = reader1->GetOutput();
    auto rb2v = reader2->GetOutput();
    auto rb3v = reader3->GetOutput();
    auto rf_up_v = reader4->GetOutput();
    auto rf_down_v = reader5->GetOutput();

    // 依次调用每个切割函数
    // 切出尾状段
    auto liver1 = cut_tail(liver_cl, rb1v, rb2v, rb3v, rf_up_v, rf_down_v, left_inside_v, left_ou_v, left_od_v, tail,mbounds_liver);
    // 提取右半叶
    auto [liver_right,liver_left_pre] = cut_right_lobe(liver, rb1v, rb2v, rb3v, rf_up_v, rf_down_v, left_inside_v, tail,mbounds_liver);
    // 提取左半叶
    auto liver_left = cut_left_lobe(liver_left_pre, left_inside_v, left_od_v, left_ou_v, tail,mbounds_liver);
    // 右后3段
    auto liver_r1 = cut_rb3b(liver_right, rb1v, rb2v, rf_up_v, rf_down_v, rb3v, tail,mbounds_liver);
    // 右后2段
    auto liver_r2 = cut_rb2b(liver_r1, rb1v, rb3v, rf_up_v, tail, rb2v,mbounds_liver);
    // 右前上段
    auto liver_r3 = cut_rf_up(liver_r2, rb1v, rf_down_v, rf_up_v,mbounds_liver);
    // 右后1段和右前下段
    cut_rb1b_rf_down(liver_r3, rb1v, rf_down_v,mbounds_liver);
    // 左内叶段
    auto liver_l1 = cut_left_inside(liver_left, left_ou_v, left_od_v, left_inside_v,mbounds_liver);
    // 左外上&下段
    cut_left_ou_od(liver_l1, left_od_v, left_ou_v,mbounds_liver);

    return EXIT_SUCCESS;
}