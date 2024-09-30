#include "cut_rb3b.h"
#include "liver_cut.h" 
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第四刀的操作
vtkSmartPointer<vtkPolyData> cut_rb3b(vtkSmartPointer<vtkPolyData> liver_right, 
              vtkSmartPointer<vtkPolyData> rb1v, 
              vtkSmartPointer<vtkPolyData> rb2v, 
              vtkSmartPointer<vtkPolyData> rf_up_v, 
              vtkSmartPointer<vtkPolyData> rf_down_v, 
              vtkSmartPointer<vtkPolyData> rb3v, 
              vtkSmartPointer<vtkPolyData> tail,
              double (&mbounds)[6]) 
{
    // 第4刀，右后3段
    auto exc_rb3v = mergePolyData(rb1v, rb2v, rf_down_v, rf_up_v, tail);
    auto [rb3b, rb3b_points] = find_points(exc_rb3v, rb3v);
    auto imp4 = plane_generator(rb3b);
    auto [liver_r1, rb3b_poly] = cut_run2(imp4, liver_right,mbounds);

    // 保存右后3段切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rb3b.stl");
    stlWriter->SetInputData(rb3b_poly);
    stlWriter->Write();

    vtkSmartPointer<vtkSampleFunction> sampleFunction = vtkSmartPointer<vtkSampleFunction>::New();
    sampleFunction->SetImplicitFunction(imp4); // 您的隐式函数
    sampleFunction->SetModelBounds(mbounds[0], mbounds[1], mbounds[2], mbounds[3], mbounds[4], mbounds[5]); // 设置采样范围
    sampleFunction->SetSampleDimensions(50, 50, 50); // 设置采样分辨率
    sampleFunction->ComputeNormalsOff();
    sampleFunction->Update();

    vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->SetInputConnection(sampleFunction->GetOutputPort());
    contourFilter->SetValue(0, 0.0); // 提取等值为 0 的面，即隐式函数的零等值面
    contourFilter->Update();

    stlWriter->SetFileName("C:/code/liver_cut/new_imp/imp_rb3b.stl"); // 设置输出文件名
    stlWriter->SetInputConnection(contourFilter->GetOutputPort());
    stlWriter->Write();


    std::cout << "Successfully saved rb3b cut" << std::endl;
    std::cout << "rb3b cut completed" << std::endl;

    return liver_r1;
}
