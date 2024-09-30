#include "cut_left_lobe.h"
#include "liver_cut.h"  // 包含你已经模块化的函数，比如 plane_generator, cut_run2 等
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第三刀的操作
vtkSmartPointer<vtkPolyData> cut_left_lobe(vtkSmartPointer<vtkPolyData> liver_left_pre, 
                   vtkSmartPointer<vtkPolyData> left_inside_v, 
                   vtkSmartPointer<vtkPolyData> left_od_v, 
                   vtkSmartPointer<vtkPolyData> left_ou_v, 
                   vtkSmartPointer<vtkPolyData> tail,
                    double (&mbounds)[6]) 
{
    // 第3刀，提取左半叶
    auto left_side_v = mergePolyData(left_inside_v, left_od_v, left_ou_v);
    auto [cut_left_side, cut_left_side_points] = find_points(left_side_v, tail);
    auto imp3 = plane_generator(cut_left_side);
    auto [liver_left, t] = cut_run2(imp3, liver_left_pre,mbounds);

    // 保存左半叶切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/liver_left.stl");
    stlWriter->SetInputData(liver_left);
    stlWriter->Write();

    vtkSmartPointer<vtkSampleFunction> sampleFunction = vtkSmartPointer<vtkSampleFunction>::New();
    sampleFunction->SetImplicitFunction(imp3); // 您的隐式函数
    sampleFunction->SetModelBounds(mbounds[0], mbounds[1], mbounds[2], mbounds[3], mbounds[4], mbounds[5]); // 设置采样范围
    sampleFunction->SetSampleDimensions(50, 50, 50); // 设置采样分辨率
    sampleFunction->ComputeNormalsOff();
    sampleFunction->Update();

    vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->SetInputConnection(sampleFunction->GetOutputPort());
    contourFilter->SetValue(0, 0.0); // 提取等值为 0 的面，即隐式函数的零等值面
    contourFilter->Update();

    stlWriter->SetFileName("C:/code/liver_cut/new_imp/imp_left_lobe.stl"); // 设置输出文件名
    stlWriter->SetInputConnection(contourFilter->GetOutputPort());
    stlWriter->Write();


    std::cout << "Successfully saved liver_left cut" << std::endl;
    std::cout << "Liver_left cut completed" << std::endl;

    return liver_left;
}
