#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPolyData.h>
#include <Eigen/Dense>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <string>

// 计算单个 vtkPolyData 的质心
Eigen::Vector3d cal_centroid(const vtkSmartPointer<vtkPolyData>& segmentPolyData) {
    Eigen::Vector3d overallCentroid(0.0, 0.0, 0.0);

    // 检查 segmentPolyData 是否有效
    if (!segmentPolyData || !segmentPolyData->GetPoints()) {
        // 返回零向量或根据需要处理错误
        return overallCentroid;
    }

    vtkIdType totalPoints = segmentPolyData->GetNumberOfPoints();

    if (totalPoints > 0) {
        for (vtkIdType i = 0; i < totalPoints; ++i) {
            double p[3];
            segmentPolyData->GetPoint(i, p);
            overallCentroid += Eigen::Vector3d(p[0], p[1], p[2]);
        }
        overallCentroid /= static_cast<double>(totalPoints);
    } else {
        // 处理没有点的情况
        // 返回零向量或根据需要处理错误
    }

    return overallCentroid;
}

// 实现炸裂效果的函数
void blow_up(std::unordered_map<std::string, vtkSmartPointer<vtkPolyData>>& segmentsPolyData, const Eigen::Vector3d& liver_center, double scaleFactor) {
    for (auto& pair : segmentsPolyData) {
        const std::string& segmentName = pair.first;
        vtkSmartPointer<vtkPolyData>& segmentPolyData = pair.second;

        vtkPoints* points = segmentPolyData->GetPoints();
        if (!points) {
            continue; // 跳过空指针
        }

        vtkIdType numPoints = points->GetNumberOfPoints();
        if (numPoints == 0) {
            continue; // 跳过没有点的分段
        }

        // 计算分段质心
        Eigen::Vector3d segmentCentroid(0.0, 0.0, 0.0);
        for (vtkIdType i = 0; i < numPoints; ++i) {
            double p[3];
            points->GetPoint(i, p);
            segmentCentroid += Eigen::Vector3d(p[0], p[1], p[2]);
        }
        segmentCentroid /= static_cast<double>(numPoints);

        // 计算位移向量并缩放
        Eigen::Vector3d moveVector = (segmentCentroid - liver_center) * scaleFactor;

        // 直接在 vtkPoints 中修改点的位置
        for (vtkIdType i = 0; i < numPoints; ++i) {
            double p[3];
            points->GetPoint(i, p);

            // 使用 Eigen 进行向量运算
            Eigen::Vector3d newPoint = Eigen::Vector3d(p[0], p[1], p[2]) + moveVector;

            // 更新点的位置
            points->SetPoint(i, newPoint[0], newPoint[1], newPoint[2]);
        }

        // 通知 VTK 数据已更新
        segmentPolyData->Modified();
    }
}

int main() {
    // 创建 STL 阅读器
    // 为肝脏整体模型创建单独的阅读器
    vtkSmartPointer<vtkSTLReader> liverReader = vtkSmartPointer<vtkSTLReader>::New();
    liverReader->SetFileName("C:/code/liver_cut/new_stl/liver.stl");
    liverReader->Update();
    vtkSmartPointer<vtkPolyData> liver = vtkSmartPointer<vtkPolyData>::New();
    liver->DeepCopy(liverReader->GetOutput());

    // 创建分段的哈希表（使用 unordered_map）
    std::unordered_map<std::string, vtkSmartPointer<vtkPolyData>> segmentsPolyData;

    // 定义分段名称和对应的文件名
    std::vector<std::pair<std::string, std::string>> segmentFiles = {
        {"tail", "C:/code/liver_cut/new_stl/tail.stl"},
        {"rb1b", "C:/code/liver_cut/new_stl/rb1b.stl"},
        {"rb2b", "C:/code/liver_cut/new_stl/rb2b.stl"},
        {"rb3b", "C:/code/liver_cut/new_stl/rb3b.stl"},
        {"rf_up", "C:/code/liver_cut/new_stl/rf_up.stl"},
        {"rf_down", "C:/code/liver_cut/new_stl/rf_down.stl"},
        {"left_inside", "C:/code/liver_cut/new_stl/left_inside.stl"},
        {"left_ou", "C:/code/liver_cut/new_stl/left_ou.stl"},
        {"left_od", "C:/code/liver_cut/new_stl/left_od.stl"}
    };

    // 读取每个分段的 STL 文件并存入哈希表
    for (const auto& pair : segmentFiles) {
        const std::string& segmentName = pair.first;
        const std::string& fileName = pair.second;

        vtkSmartPointer<vtkSTLReader> segmentReader = vtkSmartPointer<vtkSTLReader>::New(); // 为每个文件创建新的阅读器
        segmentReader->SetFileName(fileName.c_str());
        segmentReader->Update();
        auto segmentPolyData = vtkSmartPointer<vtkPolyData>::New();
        segmentPolyData->DeepCopy(segmentReader->GetOutput()); // 深拷贝数据

        segmentsPolyData[segmentName] = segmentPolyData;
    }

    // 计算肝脏整体的质心
    Eigen::Vector3d liver_center = cal_centroid(liver);

    // 设置缩放因子，可以根据需要调整
    double scaleFactor = 1.0;

    // 应用炸裂效果
    blow_up(segmentsPolyData, liver_center, scaleFactor);

    // 保存修改后的分段到新的路径
    for (const auto& pair : segmentsPolyData) {
        const std::string& segmentName = pair.first;
        vtkSmartPointer<vtkPolyData> segmentPolyData = pair.second;

        vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
        std::string filename = "C:/code/liver_cut/blew_up/" + segmentName + ".stl";
        writer->SetFileName(filename.c_str());
        writer->SetInputData(segmentPolyData);
        writer->Write();
    }

    return EXIT_SUCCESS;
}
