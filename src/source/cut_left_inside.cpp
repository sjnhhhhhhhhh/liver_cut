#include "cut_left_inside.h"
#include "liver_cut.h"  // 包含你已经模块化的函数，比如 plane_generator, cut_run2 等
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第八刀的操作
vtkSmartPointer<vtkPolyData> cut_left_inside(vtkSmartPointer<vtkPolyData> liver_left, 
                     vtkSmartPointer<vtkPolyData> left_ou_v, 
                     vtkSmartPointer<vtkPolyData> left_od_v, 
                     vtkSmartPointer<vtkPolyData> left_inside_v,
                     double (&mbounds)[6]) 
{
    // 第8刀，左内叶段
    auto exc_left_inside_v1 = mergePolyData(left_ou_v, left_od_v);
    auto [cut_left_inside1, cut_left_inside_points1] = find_points(exc_left_inside_v1, left_inside_v);
    auto imp8 = plane_generator(cut_left_inside1);
    auto [left_inside_poly, liver_l1] = cut_run2(imp8, liver_left,mbounds);

    // 保存左内叶段切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/left_inside.stl");
    stlWriter->SetInputData(left_inside_poly);
    stlWriter->Write();

    vtkSmartPointer<vtkSampleFunction> sampleFunction = vtkSmartPointer<vtkSampleFunction>::New();
    sampleFunction->SetImplicitFunction(imp8); // 您的隐式函数
    sampleFunction->SetModelBounds(mbounds[0], mbounds[1], mbounds[2], mbounds[3], mbounds[4], mbounds[5]); // 设置采样范围
    sampleFunction->SetSampleDimensions(50, 50, 50); // 设置采样分辨率
    sampleFunction->ComputeNormalsOff();
    sampleFunction->Update();

    vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->SetInputConnection(sampleFunction->GetOutputPort());
    contourFilter->SetValue(0, 0.0); // 提取等值为 0 的面，即隐式函数的零等值面
    contourFilter->Update();

    stlWriter->SetFileName("C:/code/liver_cut/new_imp/imp_left_inside.stl"); // 设置输出文件名
    stlWriter->SetInputConnection(contourFilter->GetOutputPort());
    stlWriter->Write();


    std::cout << "Successfully saved left_inside cut" << std::endl;
    std::cout << "left_inside cut completed" << std::endl;

    return liver_l1;
}
