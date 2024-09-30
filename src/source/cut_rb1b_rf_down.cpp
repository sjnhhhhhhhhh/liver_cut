#include "cut_rb1b_rf_down.h"
#include "liver_cut.h"  // 包含你已经模块化的函数，比如 plane_generator, cut_run2 等
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第七刀的操作
void cut_rb1b_rf_down(vtkSmartPointer<vtkPolyData> liver_r3, 
                      vtkSmartPointer<vtkPolyData> rb1v, 
                      vtkSmartPointer<vtkPolyData> rf_down_v,
              double (&mbounds)[6]) 
{
    // 第7刀，右后1段与右前下段
    auto [rb1b, rb1b_points] = find_points(rf_down_v, rb1v);
    auto imp7 = plane_generator(rb1b);
    auto [rf_down_poly, rb1b_poly] = cut_run2(imp7, liver_r3,mbounds);

    // 保存右后1段切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rb1b.stl");
    stlWriter->SetInputData(rb1b_poly);
    stlWriter->Write();
    std::cout << "rb1b cut completed" << std::endl;

    vtkSmartPointer<vtkSampleFunction> sampleFunction = vtkSmartPointer<vtkSampleFunction>::New();
    sampleFunction->SetImplicitFunction(imp7); // 您的隐式函数
    sampleFunction->SetModelBounds(mbounds[0], mbounds[1], mbounds[2], mbounds[3], mbounds[4], mbounds[5]); // 设置采样范围
    sampleFunction->SetSampleDimensions(50, 50, 50); // 设置采样分辨率
    sampleFunction->ComputeNormalsOff();
    sampleFunction->Update();

    vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->SetInputConnection(sampleFunction->GetOutputPort());
    contourFilter->SetValue(0, 0.0); // 提取等值为 0 的面，即隐式函数的零等值面
    contourFilter->Update();

    stlWriter->SetFileName("C:/code/liver_cut/new_imp/imp_rb1b.stl"); // 设置输出文件名
    stlWriter->SetInputConnection(contourFilter->GetOutputPort());
    stlWriter->Write();


    // 保存右前下段切割后的 STL 文件
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rf_down.stl");
    stlWriter->SetInputData(rf_down_poly);
    stlWriter->Write();
    std::cout << "rf_down cut completed" << std::endl;

}
